- 下载源码，ndk-build可执行文件，更改ceserver名为hao（防止检测）。https://blog.csdn.net/qq_43355637/article/details/126941992
- 开启：
  ```
  adb shell
  su
  cd /data/local/tmp
  ./hao
  ```
  转发端口：
  ```
  adb forward tcp:52736 tcp:52736
  ```
