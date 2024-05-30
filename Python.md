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

