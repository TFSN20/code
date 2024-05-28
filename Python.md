- 安装Python解释器时，勾选添加环境变量会将D:\ProgramFiles\Codes\Python\Python312\和D:\ProgramFiles\Codes\Python\Python312\Scripts路径添加到环境变量
- 版本 `python3 -V`（区分大小）或 `python -V`
- 查看python安装目录
```
python
import sys
print(sys.prefix)
```
- 查看pip版本 `pip —versaion` 和 `pip3 --version` 一样
- 查看包 ``pip show xxx``
- 查看所有第三方包名称 `pip list`
- 查看第三方包位置，一般位于python安装目录下的Lib文件夹下的site-packages文件夹下
- 查看pip镜像源
```
pip config list -v
```
- 查看标准库路径 `python` ，一般是python安装目录下的Lib文件夹

```
python
import os 
print(os.**file**)
```
- 查看引用所有引用模块路径
```
python
import sys
print(sys.path)
```
