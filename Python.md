- 安装Python解释器时，勾选添加环境变量会将D:\ProgramFiles\Codes\Python\Python312\和D:\ProgramFiles\Codes\Python\Python312\Scripts路径添加到环境变量
- 修改全局pip源，更改C:\ProgramData\pip\pip.ini(没有则创建)内容为（注意网址http后面有s，否则下载不了）
  ```
  [global]
  index-url = https://mirrors.aliyun.com/pypi/simple/
  ```
- 版本查看 `python3 -V`（区分大小）或 `python -V`或`python --version`
- 查看python安装目录
  cmd键入python进入解释器模式
  ```
  import sys
  print(sys.prefix)
  ```
- 查看pip版本和文件夹位置 `pip —versaion` 或 `pip3 --version` 或`pip -V`
- 查看第三方包版本和本地位置 ``pip show xxx``
- 查看所有第三方包名称和版本 `pip list`
- 查看第三方包位置，一般位于python安装目录下的Lib文件夹下的site-packages文件夹下，这也是pip的位置所在，所以pip也是一个第三方包
- 查看pip镜像源（只是列举配置文件的本地位置，不一定存在）
  ```
  pip config list -v
  ```
- 标准库路径一般是python安装目录下的Lib文件夹
- 查看引用所有引用模块路径（里面就有标准库路径）
  ```
  import sys
  print(sys.path)
  ```
- 安装特定依赖文件：如果你有一个requirements.txt文件，存有所有需要的包及其版本，可以使用以下命令进行安装：`pip install -r requirements.txt`
- 在Windows上可以通过更改环境变量的顺序来选择需要的Python版本
  ```
  D:\ProgramFiles\Codes\Python\Python312\Scripts\
  D:\Programfiles\Codes\Python\Python312\
  ```
  Scripts用于pip版本，这个可以保留一个一直提前，这样就不用切换python版本时再下载包，对包的操作还是切换前pip对包的操作
## pytorch
- 安装:
  根据[官网](https://pytorch.org/)python版本和是否有英伟达显卡（没有选择则选择cpu）复制pip下载链接。
- tensorboard
  - 安装```pip install tensorboard```
  - 导入SummaryWriter显示：未存取“SummaryWriter”Pylance未从模块“torch.utils.tensorboard”导出“SummaryWriter”
    改为从"torch.utils.tensorboard.writer"导入。
    ```
    from torch.utils.tensorboard import SummaryWriter  # type: ignore
    ```
    不用管，可以正常运行
  - vsc终端或cmd指定目录```tensorboard --logdir C:\Users\Administrator\Desktop\1```
  - 指定端口```tensorboard --logdir=C:\Users\Administrator\Desktop\1 --port=6007```
  - 关闭端口（即关闭tensorboard）：终端键入ctrl+c中断
  - cv2读取显示
    ```
    from torch.utils.tensorboard import SummaryWriter  # type: ignore
    from torchvision import transforms
    import cv2
    
    img = cv2.imread(r"D:\Downloads\Compressed\hymenoptera_data\hymenoptera_data\train\ants\0013035.jpg") #默认3通道读取
    img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB) #cv2图像数据默认是BGR
    tensor = transforms.ToTensor()
    img_tensor = tensor(img)
    norm = transforms.Normalize([0.5, 0.5, 0.5], [0.5, 0.5, 0.5])
    img_norm = norm(img_tensor)
    writer = SummaryWriter(r'C:\Users\Administrator\Desktop\1')
    writer.add_image('norm', img_tensor, 1) # 使用tensor类型作为参数时，shape和dataformats默认参数一样都是'CHW'
    writer.add_image('norm', img_norm, 2) # 2作为step参数，可以是任意数字，显示时会按照数值大小（取整-2.6取-2）排列显示.
    writer.close()
    ```
  - transforms.Normalize这一步通常在预处理阶段进行，目的是使模型训练更加稳定和快速。
  - cv2图像shape为高宽通道，tensor类型为通道高宽，SummaryWriter().add_image()使用tensor类型图片时，不用设置默认的dataformats的参数'CHW'，故推荐。
    ```
    import numpy as np
    from PIL import Image

    image = Image.open('path/to/your/image.jpg')
    image_array = np.array(image)
    print(image_array.shape)
    ```
    PIL图像数据为自己独有的类型，SummaryWriter().add_image()只能使用numpy和tensor数据类型，转为numpy类型后shape为高宽通道，还需要增加dataformats='HWC'参数，好在PIL是自带的库。
  - Resize：transforms.Resize((100))表示最短边缩放到100，长边等比缩放，但是100不能大于最短边边长；transforms.Resize((800,800))表示图像高宽缩放到此像素，没有限制。
  - RandomCrop：一个参数时对应的不能大于对应边的像素，一个参数情况包括512和(512,512)【高宽】两种情况，两个参数时可以大于，且即使transforms.RandomCrop(1, 1)也会取到图片外的黑色像素点。
## os模块
- 获取目录下的子文件的文件全名，以数组形式
  ```
  os.listdir(r'd:\Downloads\Compressed\hymenoptera_data\hymenoptera_data\train')
  # 得到数组[test.html, 1.mp4, stamp]
  ```
- 连接路径
  ```
  os.path.join('D:/','a','b')
  ```
## 常用
- tensor_img[:, 0, 0]
假设 tensor_img 是一个形状为 (3, H, W) 的Tensor，其中 3 是通道数（分别对应于RGB），H 是图像的高度，W 是图像的宽度。
tensor_img 代表整个图像的Tensor。
tensor_img[:, 0, 0] 是对这个Tensor进行索引操作，具体解释如下：
: 表示选择所有的通道。
0 表示高度的第一个像素（从上到下的第一个像素）。
0 表示宽度的第一个像素（从左到右的第一个像素）。

