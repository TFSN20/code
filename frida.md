- 注意python, frida, frida-tools, frida-server的版本和android架构
  ```
  adb shell getprop ro.product.cpu.abi # 架构查看
  ```
  - 下载frida-server：https://github.com/frida/frida/releases
  - 一个对应：frida-16.4.8，frida-tools-12.5.0
- 查看tcp为5037的占用情况
  ```
  netstat -ano | findstr :5037
  ```
  根据PID杀死进程
  ```
  taskkill /PID 进程的PID /F
  ```
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
  - 查看权限：
  ```
  adb shell ls -l /data/local/tmp/frida-server
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
## 常用脚本
- 开启frida：
  ```
  adb shell
  su
  cd /data/local/tmp
  ./av -l 127.0.0.1:27777 &
  ```
- 查看android frida端口
  ```
  adb kill-server
  adb shell ps | findstr av
  ```
  杀死该端口：
  ```
  adb shell su -c "kill -9 9820"
  ```
- 转发端口：
  ```
  adb forward tcp:27777 tcp:27777
  netstat -aon | findstr :27777
  frida-ps -U
  ```
- js函数
  ```
  let fn = function () {
      console.log('fn will java perform')
      Java.perform(function () {
          console.log('fn will hook');
  
          var WebView = Java.use("android.webkit.WebView");
  
          WebView.loadUrl.overload("java.lang.String").implementation = function (s) {
  
              this.loadUrl.overload("java.lang.String").call(this, s);
          };
  
      })
  }
  ```
- python启动
  ```
  import frida
  import sys
  
  PACKAGE = 'com.greenpoint.android.mc10086.activity'
  HOST = '127.0.0.1:27777'
  
  def on_message(message,data):
      if message['type'] == 'send':
          print(message['payload'])
      elif message['type'] == 'error':
          print(message['stack'])
  
  if __name__ == '__main__':
      try:
  
          manager = frida.get_device_manager()
          device = manager.add_remote_device(HOST)
          jscode = open(r'D:\cal\1.js','r',encoding='UTF-8').read()
          pid = device.spawn([PACKAGE])
          session = device.attach(pid)
  
          '''
          attach模式
          '''
          # manager = frida.get_device_manager()
          # device= manager.add_remote_device(HOST)
          # jscode = open(r'D:\cal\1.js','r',encoding='UTF-8').read()
          # session = device.attach('中国移动')
          
  
          script = session.create_script(jscode)
          script.on('message',on_message)
          script.load()

        script.exports_sync.hook_webvi()

        device.resume(pid)

        sys.stdin.read()
    except Exception as e:
        print(e)
        print('frida hook error')
  ```
- 查看app pid（运行状态）
  ```
  adb shell ps | findstr package_name
  ```
- 查看app进程ID（PID）
  ```
  adb shell
  su
  ls /proc/pid/task
  ```
- 查看maps 过滤
  ```
  cat /proc/pid/maps | grep frida
  ```
- 线程ID（TID）：
  ```
  cat /proc/pid/task/tid/status
  ```
- tid全部的name status:
  ```
  for i in $(ls /proc/11157/task/); do
    echo "$i"
    head -1 /proc/11157/task/$i/status
    echo "—"
  done
  ```
- 终端启动（spawn）
  ```
  frida -l "D:\OneDrive\Codes\frida\1\frida_script.js" -H 127.0.0.1:27777 -o "C:\Users\Administrator\Desktop\1.txt" -f io.cordova.zhqy
  ```
  ```
  frida -U -f tv.danmaku.bili -l "D:\OneDrive\Codes\frida\1\frida_script.js"
  ```

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
- 有些app会向firda-server发送请求，如果回应为REJECT字符，则检测到，使用strcmp或strstr比较两个字符指针里的字符（回应字符和REJECT字符）是否一样，一样则返回0，所以可以hook这两个函数，如果函数的前两个参数中有一个是REJECT字符，证明app有此项检测，我们只需返回非0即可。（先hook再运行app可以保证即使app对strcmp进行了hook，也可被修改hook。app是否可以在运行时通过网络请求js，js提供比较函数使用，或则让frida没法hook 其他平台？接口提供的比较函数）
  ```
  strcmp("A", "B") // <0
  strcmp("REJEC", "REJECT") // 亦<0，所以不能
  strcmp("AB", "B") // <0 A小于B的ASCII码, 只比第一位，一样则往后比较，没有后一位则\0为终止符比较，ASCII码最小为0，即\0
  ```
  ```
  strstr('REJECT', 'REJECT') // 返回指针
  strstr('REJECT', 'REJECT1') // NULL frida js 使用0
  strstr('REJECT', 'JECT') // 返回指针
  ```
- 文件描述符（file descriptors）检测，每个文件描述符是一个符号链接，指向它所引用的文件或资源：frida-server运行时会涉及创建和管理符号链接，"/proc/self/fd/"目录下记录了这些链接，当链接中包含可疑字符时即检测到。这些链接以\开头，所以只需要strstr第二个参数（可疑字符）判断是否在第一个参数（链接）内即可。
- map文件检测：`/proc/self/maps` 是一个特殊的文件，它包含了当前进程的内存映射信息。当你打开这个文件时，它会显示一个列表，其中包含了进程中每个内存区域的详细信息。其中包含frida-xxx.so等文件（文件路径（如果该内存区域映射了一个文件））。
## 无法调试和一些问题
- frida调试时，有些app会停在启动页面，这时候的hook代码都是不起作用的！！！（有些还会重启手机）这时候就需要使用终端模式注入了，速度更快，因为一些app整体代码比较简单，加载反调试较为靠前。
- 调用rpc函数时
  ```
  result = script.exports_sync.greet('world')
  aaa = script.exports_sync.asyncfetch()
  print(f"Greet result: {result}")
  print(f"async result: {aaa}")
  ```
  正确，但是
  ```
  print(f"Greet result: {script.exports_sync.greet('world')}")
  print(f"async result: {script.exports_sync.asyncfetch()}")
  ```
  或者
  ```
  print(f"async result: {script.exports_sync.asyncfetch()}")
  print(f"Greet result: {script.exports_sync.greet('world')}")
  ```  
  都错误：script has been destroyed，因为调用一次以上的函数时，异步函数在python的模板字符串中不行。
  
