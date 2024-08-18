- 下载源码，ndk-build可执行文件，或则https://www.cheatengine.org/downloads.php直接下载，更改ceserver名为hao（防止检测）。https://blog.csdn.net/qq_43355637/article/details/126941992
- 开启：
  ```
  adb shell
  su
  cd /data/local/tmp
  ./hao -p 28888
  ```
  新开终端，转发端口：
  ```
  adb forward tcp:28888 tcp:28888
  netstat -aon | findstr :28888
  ```
  防止检测：edit->settings->Debugger Options->Try to prevent detection of the debugger
  
  ce：file->open process->network host:127.0.0.1 port:28888
- 字节（Byte）是计算机中存储和处理数据的基本单位。一个字节等于8位（bits），每位可以是0或1。字节是计算机内存中存储数据的基本构建块，用来表示数据的基本单位。

  以下是一些常见数据类型及其在内存中的字节大小：
  8位（1字节）：通常用于存储一个字符（如ASCII码），或是一个小的整数（0到255）。
  
  16位（2字节）：用于存储较小的整数（-32768到32767），也用于一些短字符编码。
  
  32位（4字节）：常用于存储较大的整数（-2147483648到2147483647），或是单精度浮点数（float）。
  
  64位（8字节）：用于存储更大的整数（-9223372036854775808到9223372036854775807），或是双精度浮点数（double）。
  
  在游戏中存储坐标的字节数
  
  单精度浮点数（float）：通常占用 4字节。这是游戏中存储坐标最常见的数据类型，因为它能提供足够的精度，同时又不会占用过多的内存。
  
  双精度浮点数（double）：占用 8字节。用于需要更高精度的场景。
