## 颜色相似度
要计算两种颜色之间的径向相似度（颜色深浅）和圆周相似度（颜色品种），首先需要明白径向相似度指的是颜色的亮度或深度的相似度，而圆周相似度指的是颜色本身的色相或品种的相似度。这意味着我们需要同时比较两种颜色的亮度和色相。
在进行计算之前，应该将RGB颜色值转换为一个能分离色相、饱和度和亮度的颜色空间，如HSL或HSV，因为它们更适合进行这样的比较。这里我们以HSL颜色空间为例进行解释。
RGB转HSL：
使用颜色colorA=[123,135,212]和colorB=[12,34,23]，首先需要将它们转换成HSL颜色空间。这个转换可能需要复杂的公式，但许多编程语言都有内置函数或库来完成这个工作。
圆周相似度（色相）计算：
比较每种颜色的色相成分（色相为色轮上的角度）。这可以通过找出两个色相角度的差异，并将其规范到区间[0, 180]，因为色轮上色相的最大差异是180度。
python
Copy
def hue_difference(hue1, hue2):
    diff = abs(hue1 - hue2)
    return min(diff, 360 - diff) # 规范到区间[0,180]
径向相似度（饱和度/亮度）计算：
对于径向相似度，比较两种颜色的饱和度和亮度值。通过互相减去亮度值及饱和度值并取其绝对值可以得出亮度和饱和度的相似度。
python
Copy
def radial_difference(lightness1, lightness2, saturation1, saturation2):
    lightness_diff = abs(lightness1 - lightness2)
    saturation_diff = abs(saturation1 - saturation2)
    return lightness_diff, saturation_diff
结合相似度：
根据应用的不同，你可能需要不同的方式来权衡色相差异和亮度/饱和度差异。为了得出综合的相似度得分，你可以使用加权求和方法，其中你可以根据色相和亮度/饱和度的重要性来分配权重。
以下是一个结合了所有这些步骤的简单Python脚本示例：
```
import colorsys

def rgb_to_hsl(r, g, b):
    return colorsys.rgb_to_hls(r/255.0, g/255.0, b/255.0)

def calc_similarity(colorA, colorB):
    hslA = rgb_to_hsl(*colorA)
    hslB = rgb_to_hsl(*colorB)

    hue_diff = hue_difference(hslA[0]*360, hslB[0]*360)
    lightness_diff, saturation_diff = radial_difference(hslA[1], hslB[1], hslA[2], hslB[2])

    # 将色相、饱和度和亮度差异合并成一个相似度指标（示例）
    combined_similarity = (1 - hue_diff/180) * (1 - saturation_diff) * (1 - lightness_diff)

    return {
        'hue_similarity': 1 - hue_diff/180,
        'saturation_similarity': 1 - saturation_diff,
        'lightness_similarity': 1 - lightness_diff,
        'combined_similarity': combined_similarity
    }
```
# 使用所提供的RGB颜色示例
```
colorA = [123,135,212]
colorB = [12,34,23]

similarity = calc_similarity(colorA, colorB)
print(similarity)
```
请注意，这个脚本假设了一种简化的相似度组合方式，在你的应用中，你可能需要调整这些相似度的组合方式。在色彩科学中，存在更复杂的模型和距离度量，例如CIEDE2000，它们更准确地量化了人类视觉所感知的颜色差异。
