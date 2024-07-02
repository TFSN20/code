在PyTorch中,当使用torch.cat()函数进行张量拼接时, dim参数确定了沿着哪个维度进行拼接,以下是
参数的详细解释：
1. dim=0 （第一个维度）：
会在每个张量的第一个维度上进行拼接，
对于3D张量而言，这个维度指的是张量的“块”数量。
示例:如果有两个3×3×4形状的张量,通过在 dim-e 拼接会得到一个6×3×4的张量,x
2. dim=1 （第二个维度）：
会在每个张量的第二个维度上进行拼接。
对于3D张量而言，这个维度指的是每个“块”中的行数。
示例:如果有两个3x3×4形状的张量,通过在dim-1 拼接会得到一个3×6×4的张量,x
3. dim=-1 或 dim=2 (第三个维度):
会在每个张量的第三个维度上进行拼接,
对于3D张量而言，这个维度指的是每一行中的列数。
示例：如果有两个3× 3× 4形状的张量，通过在 dim--1 或 dim-2 拼接会得到一个3 ×3× 8 的张量。