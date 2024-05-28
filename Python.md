- 安装Python解释器时，勾选添加环境变量会将D:\ProgramFiles\Codes\Python\Python312\和D:\ProgramFiles\Codes\Python\Python312\Scripts路径添加到环境变量
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
- 查看pip镜像源
  ```
  pip config list -v
  ```
- 标准库路径一般是python安装目录下的Lib文件夹
- 查看引用所有引用模块路径（里面就有标准库路径）
  ```
  import sys
  print(sys.path)
  ```
