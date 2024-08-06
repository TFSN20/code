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
  - 启动：
  ```
  cd /data/local/tmp/
  ./frida-server &
  一般指定端口启动：./frida-server -l 127.0.0.1:27777 &
  adb forward tcp:27777 tcp:27777 # 端口转发
  frida-ps -U # 另起cmd输入
  ```
  - 查看端口占用
  ```
  netstat -aon | findstr :27042
  ```
  - 根据进程PID杀死
  ```
  tasklist /FI "PID eq 12300"
  ```
  - 一键运行：
  ```
  adb shell
  su
  cd /data/local/tmp
  ./frida-server -l 127.0.0.1:27777 &
  ```
  - 一键转发
  ```
  adb forward tcp:27777 tcp:27777
  netstat -aon | findstr :27777
  frida-ps -U
  ```
- 两种方式hook
  - spawn模式：重新启动hook，hook的时机非常早（即在App启动阶段）
    ```
    device = frida.get_remote_device()
    pid = device.spawn(['com.icbc'])
    device.resume(pid)
    time.sleep(1)
    session = device.attach(pid)
    ```
  - attach模式：附加hook，需要App处于启动状态
    ```
    session = frida.get_remote_device().attach('中国工商银行')
    ```
- hook闪退：frida被检测到了
  - 改端口启动frida-server
  - frida-server改名lud
    ```
    mv frida-server lud
    ```
    ```
    adb shell
    su
    cd /data/local/tmp
    ./lud -l 127.0.0.1:27777 &
    ```
