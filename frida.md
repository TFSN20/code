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
  netstat -aon | findstr :27777
  ```
  - 根据进程PID杀死
  ```
  tasklist /FI "PID eq 12300" （好像对frida无效）
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
  - 改端口启动frida-server，frida-server改名，这里改为av
    ```
    adb shell
    su
    cd /data/local/tmp
    ./av -l 127.0.0.1:27777 &
    ```
  - frida-server去特征字段，使用https://github.com/hzzheyang/strongR-frida-android/releases下载
    ```
    adb push "d:\Downloads\Compressed\hluda-server-16.4.8-android-arm64\frida-server" /data/local/tmp/
    adb shell
    su
    chmod 777 /data/local/tmp/frida-server
    ```
  - rpc.exports函数名名字必须小写，不能有_.
## frida常见检测手段
- 打印args详细：
  ```
  var count = 0;
  // Iterate over the args to determine the count
  while (args[count] !== undefined && args[count] !== null) {
      console.log("args[" + count + "]: " + args[count].readCString());
      console.log(`args[${count}]: ${args[count].readCString()}, ${args[count]};`);
      count++;
  }
  ```
- 有些app会向firda-server发送请求，如果回应为REJECT字符，则检测到，使用strcmp或strstr比较两个字符指针里的字符（回应字符和REJECT字符）是否一样，一样则返回0，所以可以hook这两个函数，如果函数的前两个参数中有一个是REJECT字符，证明app有此项检测，我们只需返回非0即可。（但是否有这样的情况，如果strcmp("REJEC", "REJECT")）
  ```
  if(strlen(args[0].readCString())==6){
    if(strcmp("REJEC", "REJECT"))
  }
  ```
  ```
  strcmp("A", "B") // <0
  strcmp("REJEC", "REJECT") // 亦<0，所以不能
  strcmp("AB", "B") // <0 A小于B的ASCII码, 只比第一位，一样则往后比较，没有后一位则\0为终止符比较，ASCII码最小为0，即\0
- 文件描述符（file descriptors）检测，每个文件描述符是一个符号链接，指向它所引用的文件或资源：frida-server运行时会涉及创建和管理符号链接，"/proc/self/fd/"目录下记录了这些链接，当链接中包含可疑字符时即检测到。这些链接以\开头，所以只需要strstr第二个参数（可疑字符）判断是否在第一个参数（链接）内即可。
  
