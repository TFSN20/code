- 注意python, frida, frida-tools, frida-server的版本和android架构
  ```
  adb shell getprop ro.product.cpu.abi # 架构查看
  ```
  - 下载frida-server：https://github.com/frida/frida/releases
  - 一个对应：frida-16.4.8，frida-tools-12.5.0
- 库安装：
  ```
  pip install frida
  pip install frida-tools
  ```
  - 查看frida库版本
  ```
  frida --version
  ```
- 推送到手机端（解压文件改名为frida-server）：
  ```
  adb push "D:\Downloads\frida-server-16.4.8-android-arm64\frida-server" /data/local/tmp/
  ```
  - 最高权限：
  ```
  adb shell
  su
  chmod 777 /data/local/tmp/frida-server
  ```
  - 测试：
  ```
  cd /data/local/tmp/
  ./frida-server
  指定端口启动：./frida-server -l 127.0.0.1:27043
  frida-ps -U # 另起cmd输入
  ```
  - 查看端口占用
  ```
  netstat -aon | findstr :10808
  ```
  - 根据进程PID杀死
  ```
  tasklist /FI "PID eq 11964"
  ```