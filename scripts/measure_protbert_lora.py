import pathlib
import urllib.request
import json
import os

def measure_protbert_lora():
    log_path = str(pathlib.Path(__file__).resolve().parent.parent / 'logs/lora_measurement.log')
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    
    url = "https://huggingface.co/Rostlab/prot_bert_bfd/raw/main/config.json"
    req = urllib.request.Request(url, headers={'User-Agent': 'AuditProtBert/1.0'})
    
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            cfg = json.loads(resp.read().decode())
            source = "HuggingFace Hub (Rostlab/prot_bert_bfd/config.json)"
    except Exception as e:
        source = f"Fallback due to {e}"
        cfg = {
            "architectures": ["BertForMaskedLM"],
            "hidden_size": 1024,
            "intermediate_size": 4096,
            "max_position_embeddings": 512,
            "num_attention_heads": 16,
            "num_hidden_layers": 30,
            "type_vocab_size": 2,
            "vocab_size": 30
        }

    hidden_size = cfg.get("hidden_size", 1024)
    num_layers = cfg.get("num_hidden_layers", 30)
    intermediate_size = cfg.get("intermediate_size", 4096)
    num_heads = cfg.get("num_attention_heads", 16)
    vocab_size = cfg.get("vocab_size", 30)
    max_pos = cfg.get("max_position_embeddings", 512)
    type_vocab = cfg.get("type_vocab_size", 2)
    
    # Exact parameter calculation
    # 1. Embeddings
    word_embed = vocab_size * hidden_size
    pos_embed = max_pos * hidden_size
    type_embed = type_vocab * hidden_size
    embed_norm = 2 * hidden_size
    total_embeddings = word_embed + pos_embed + type_embed + embed_norm
    
    # 2. Encoder Layer
    # Self-attention:
    q_proj = hidden_size * hidden_size + hidden_size
    k_proj = hidden_size * hidden_size + hidden_size
    v_proj = hidden_size * hidden_size + hidden_size
    out_proj = hidden_size * hidden_size + hidden_size
    attn_norm = 2 * hidden_size
    
    # Feed-forward (intermediate + output):
    ffn_inter = hidden_size * intermediate_size + intermediate_size
    ffn_out = intermediate_size * hidden_size + hidden_size
    ffn_norm = 2 * hidden_size
    
    layer_params = (q_proj + k_proj + v_proj + out_proj + attn_norm) + (ffn_inter + ffn_out + ffn_norm)
    total_encoder = num_layers * layer_params
    
    # 3. Pooler
    pooler = hidden_size * hidden_size + hidden_size
    
    total_base_bert = total_embeddings + total_encoder + pooler
    
    # 4. LoRA Configuration: r=8, alpha=16, targets: query, value
    r = 8
    alpha = 16
    lora_a_q = hidden_size * r
    lora_b_q = r * hidden_size
    lora_q = lora_a_q + lora_b_q
    
    lora_a_v = hidden_size * r
    lora_b_v = r * hidden_size
    lora_v = lora_a_v + lora_b_v
    
    lora_per_layer = lora_q + lora_v
    total_lora_trainable = num_layers * lora_per_layer
    
    frac_base = (total_lora_trainable / total_base_bert) * 100
    
    # If 420M nominal parameter figure is cited:
    nominal_420m = 420_000_000
    frac_nominal = (total_lora_trainable / nominal_420m) * 100
    
    log_text = f"""=== Rostlab/prot_bert_bfd LoRA Architecture Measurement Log ===
Source: {source}

Config Parameters:
  - model_type: bert
  - num_hidden_layers (L): {num_layers}
  - hidden_size (d_model): {hidden_size}
  - intermediate_size (d_ff): {intermediate_size}
  - num_attention_heads: {num_heads}
  - vocab_size: {vocab_size}
  - max_position_embeddings: {max_pos}
  - type_vocab_size: {type_vocab}

Exact Base Parameter Breakdown:
  - Embeddings: {total_embeddings:,}
  - Per Encoder Layer: {layer_params:,}
  - Total Encoder ({num_layers} layers): {total_encoder:,}
  - Pooler Layer: {pooler:,}
  - Total Base BertModel Parameters: {total_base_bert:,} (~{total_base_bert/1e6:.2f}M)

LoRA Adaptation (r={r}, alpha={alpha}, target_modules=['query', 'value']):
  - Adapters per layer:
      * query (A: {hidden_size}x{r}, B: {r}x{hidden_size}) = {lora_q:,} params
      * value (A: {hidden_size}x{r}, B: {r}x{hidden_size}) = {lora_v:,} params
      * Total per layer: {lora_per_layer:,} params
  - Total Trainable LoRA Parameters across {num_layers} layers: {total_lora_trainable:,}
  - Trainable Fraction of Base Model ({total_base_bert:,}): {frac_base:.4f}%
  - Trainable Fraction of Nominal 420M: {frac_nominal:.4f}%

Comparison with Previous Report Assertion:
  - Previous Assertion: 1,179,648 parameters (0.28%)
  - Cause of Discrepancy: 1,179,648 corresponds to 36 layers (36 x 32,768), typical of a 36-layer architecture.
  - Measured Value: Exactly {total_lora_trainable:,} parameters across the {num_layers} layers of ProtBERT.
"""
    print(log_text)
    with open(log_path, 'w', encoding='utf-8') as f:
        f.write(log_text)
    print(f"Log saved to: {log_path}")

if __name__ == '__main__':
    measure_protbert_lora()
