# Ubuntu安装
- 不同版本有不同的安装文件夹名称，自行搜索，以下为Ubuntu22.04.3 LTS版本，D:/mklink/Ubuntu要先创建
  ```
  mklink /D "C:\Users\tfsn20\AppData\Local\Packages\CanonicalGroupLimited.Ubuntu22.04LTS_79rhkp1fndgsc" "D:/mklink/Ubuntu"
  ```
- 微软商店搜索下载（不要挂梯子），下载完成后打开（报错自行搜索错误代码，比如我的是0x8007019e，就管理员打开powershell，输入以下命令回车），新建用户名，回车，随便输入密码（界面不会有任何变化），回车。
  ```
  Enable-WindowsOptionalFeature -Online -FeatureName Microsoft-Windows-Subsystem-Linu
  ```
- 更改apt源（用的阿里源，清华源也可以，jammy是Ubuntu22.04的代号，每个大版本都不一样）
  ```
  sudo cp /etc/apt/sources.list /etc/apt/sources.list.backup #root备份源
  ```
  ```
  deb http://mirrors.aliyun.com/ubuntu/ jammy main restricted universe multiverse
  deb-src http://mirrors.aliyun.com/ubuntu/ jammy main restricted universe multiverse

  deb http://mirrors.aliyun.com/ubuntu/ jammy-updates main restricted universe multiverse
  deb-src http://mirrors.aliyun.com/ubuntu/ jammy-updates main restricted universe multiverse

  deb http://mirrors.aliyun.com/ubuntu/ jammy-backports main restricted universe multiverse
  deb-src http://mirrors.aliyun.com/ubuntu/ jammy-backports main restricted universe multiverse

  deb http://mirrors.aliyun.com/ubuntu/ jammy-security main restricted universe multiverse
  deb-src http://mirrors.aliyun.com/ubuntu/ jammy-security main restricted universe multiverse
  ```
  ```
  sudo apt update
  ```
  ```
  sudo apt upgrade
  ```
# 一些库的安装
- ```
  sudo apt-get install software-properties-common
  ```
## python库
### fenics
- 添加私人软件源
  ```
  sudo add-apt-repository ppa:fenics-packages/fenics
  ```
  下载
  ```
  sudo apt install fenics
  ```
# 常用命令
- 查看文件和文件夹
  ```
  ls
  ```
- 查看文件和文件夹的权限
  ```
  ls -la
  ```
- 创建文件
  ```
  touch filename.txt
  ```
- 删除文件夹下所有内容 先cd到目标文件夹的父目录
  ```
  rm -rf /目标文件夹
  ```
- 重命名
  ```
  mv name.py name2.py
  ```
- cd到一个文件夹a下，将所有cal开头的文件移动到文件夹a下的b文件夹
  ```
  mv cal*.* b文件夹/
  ```
- 创建一个文件
  ```
  touch xx.txt
  ```
- 编辑器编辑文本,光标移动，ctrl+x退出，y保存，回车确认文件名保存
  ```
  nano xxx.txt
  ```
- 代理全局
  ```
  export host_ip=$(ip route | grep default | awk '{print $3}')
  export http_proxy="http://${host_ip}:10808"
  export https_proxy="http://${host_ip}:10808"
  export all_proxy="socks5:/${host_ip}:10808"
  export no_proxy="localhost,127.0.0.1,*.local"
  ```
- 禁止dns自动生成 
  ```
  nano /etc/wsl.conf
  ```
  ```
  [network]
  generateResolvConf = false
  [boot]
  systemd=true
  ```
  dns修改:
  ```
  nano /etc/resolv.conf
  ```
  nameserver 223.5.5.5
  nameserver 180.76.76.76
  ```
  
  
# 常见问题
- -rw-r--r-- 1 tfsn20 tfsn20    0 Dec 16 10:32 my_script.py 能执行python3 my_script.py
  Python脚本不需要执行权限才能被运行。关键是Python解释器需要执行权限,而不是脚本本身。
