# 20231130_PCA_analysis_to_confirm_similarity_between_CCLE_HAEMATO_and_K562

# 解析の説明
K562は慢性骨髄性白血病患者から樹立された細胞株で、K562の他にも様々な白血病の患者から樹立された細胞株が存在している。
CCLEで"haematopoietic_and_lymphoid_tissue"に分類されている細胞株はK562と同じmiRNA/遺伝子発現パターンを持っているかPCA及びUMAPで解析した。

# ファイルの説明
・20231130_PCA_analysis_to_confirm_similarity_between_CCLE_HAEMATO_and_K562.R　　CCLEのmiRNA/遺伝子発現データをK562細胞株、"haematopoietic_and_lymphoid_tissue"に分類される細胞株とそれ以外の細胞株に分類し、PCA解析をおこなうことで"haematopoietic_and_lymphoid_tissue"とK562間の発現パターンの類似性を確認するスクリプト。
