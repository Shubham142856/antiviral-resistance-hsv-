import urllib.request
import json

def measure_protbert_lora():
    # ProtBERT-BFD config from HuggingFace
    url = "https://huggingface.co/Rostlab/prot_bert_bfd/raw/main/config.json"
    req = urllib.request.Request(url, headers={'User-Agent': 'AntiviralAudit/1.0 (mailto:audit@example.com)'})
    
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            cfg = json.loads(resp.read().decode())
    except Exception as e:
        print("Falling back to local / default BertConfig values:", e)
        cfg = {
            "hidden_size": 1024,
            "num_attention_heads": 16,
            "num_hidden_layers": 30,
            "intermediate_size": 4096,
            "vocab_size": 30
        }
        
    hidden_size = cfg.get("hidden_size", 1024)
    num_layers = cfg.get("num_hidden_layers", 30)
    intermediate_size = cfg.get("intermediate_size", 4096)
    vocab_size = cfg.get("vocab_size", 30)
    
    print("=== ProtBERT-BFD Architecture Config ===")
    print(f"Num hidden layers (L): {num_layers}")
    print(f"Hidden size (d_model): {hidden_size}")
    print(f"Num attention heads: {cfg.get('num_attention_heads', 16)}")
    print(f"Intermediate size: {intermediate_size}")
    print(f"Vocab size: {vocab_size}")
    
    # Calculate base transformer parameter count
    # Per layer:
    # Multi-head attention:
    # W_q: d x d + d
    # W_k: d x d + d
    # W_v: d x d + d
    # W_o: d x d + d
    # LayerNorm 1: 2 * d
    # FeedForward:
    # Linear 1: d x d_ff + d_ff
    # Linear 2: d_ff x d + d
    # LayerNorm 2: 2 * d
    
    attn_params_per_layer = 4 * (hidden_size * hidden_size + hidden_size) + 2 * hidden_size
    ffn_params_per_layer = (hidden_size * intermediate_size + intermediate_size) + (intermediate_size * hidden_size + hidden_size) + 2 * hidden_size
    layer_params = attn_params_per_layer + ffn_params_per_layer
    
    # Embeddings:
    word_embed = vocab_size * hidden_size
    pos_embed = cfg.get("max_position_embeddings", 512) * hidden_size
    token_type_embed = cfg.get("type_vocab_size", 2) * hidden_size
    embed_norm = 2 * hidden_size
    total_embed = word_embed + pos_embed + token_type_embed + embed_norm
    
    total_base_params = total_embed + (num_layers * layer_params)
    print(f"\nCalculated Total Base Model Parameters: {total_base_params:,} (~{total_base_params/1e6:.1f}M)")
    
    # LoRA parameter count for rank r=8 applied to query and value projection matrices (W_q and W_v)
    # For each adapted linear matrix of size (d_in, d_out) = (1024, 1024):
    # LoRA down-projection A: 1024 x r
    # LoRA up-projection B: r x 1024
    # Parameters per adapted matrix = 2 * 1024 * r
    
    for r in [4, 8, 16]:
        params_per_matrix = 2 * hidden_size * r
        # Adapting W_q and W_v (2 matrices per layer)
        lora_qv_per_layer = 2 * params_per_matrix
        lora_qv_total = num_layers * lora_qv_per_layer
        frac_qv = (lora_qv_total / total_base_params) * 100
        
        # Adapting W_q, W_k, W_v, W_o (all 4 attention matrices per layer)
        lora_all_attn_per_layer = 4 * params_per_matrix
        lora_all_attn_total = num_layers * lora_all_attn_per_layer
        frac_all = (lora_all_attn_total / total_base_params) * 100
        
        print(f"\n--- LoRA Rank r = {r} ---")
        print(f"Adapting W_q & W_v across {num_layers} layers:")
        print(f"  Trainable Parameters: {lora_qv_total:,}")
        print(f"  Trainable Fraction of Base Model: {frac_qv:.4f}% (~0.23% for r=8)")
        print(f"Adapting W_q, W_k, W_v, W_o across {num_layers} layers:")
        print(f"  Trainable Parameters: {lora_all_attn_total:,}")
        print(f"  Trainable Fraction of Base Model: {frac_all:.4f}%")

if __name__ == '__main__':
    measure_protbert_lora()
