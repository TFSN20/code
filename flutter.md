# 环境（绿色版本）
## JAVA
- https://adoptium.net/temurin/releases/
- Operating System: Windows  
  Architecture: x64  
  Package Type: JDK  
  Version: 17  
- 下载解压，不需要新建目录
## Android SDK
- Flutter SDK for Windows (zip)，大约2GB
- 下载解压，不需要新建目录
- 在cmdline-tools文件夹下新建latest文件夹，将android_sdk目录下所有文件移动到latest文件夹下，这是Flutter 和 Android Studio 的标准结构期望需要。
## Flutter SDK
- Android Studio 官网，下翻找到Command line tools only
- 下载zip文件，需要新建目录（如android_sdk）
## 临时环境变量配置脚本
- start_dev.bat
  ```
  @echo off
  rem =================================================================
  rem ==  Flutter 便携式开发环境启动脚本
  rem =================================================================
  
  rem 切换终端的编码为 UTF-8
  chcp 65001 > nul
  
  echo.
  echo  正在设置临时开发环境...
  echo.
  
  rem --- (可选) 设置网络代理 ---
  rem !!! 重要: 如果你的网络需要代理，请去掉下面两行的 rem 并填入你的代理地址 !!!
  set "HTTP_PROXY=http://127.0.0.1:10808"
  set "HTTPS_PROXY=http://127.0.0.1:10809"
  
  rem --- 1. 设置 Flutter SDK 的根目录 ---
  set "FLUTTER_ROOT=E:\Codes\Android\dev\flutter\flutter"
  
  rem --- 2. 设置 Android SDK 的根目录 ---
  set "ANDROID_SDK_ROOT=E:\Codes\Android\dev\android_sdk"
  
  rem --- 3. 设置 Java 开发工具包 (JDK) ---
  set "JAVA_HOME=E:\Codes\Android\dev\Temurin 17.0.17+10 - 20251031\jdk-17.0.17+10"
  
  rem --- 4. 将所有工具的 bin 目录添加到临时的 PATH 环境变量中 ---
  set "PATH=%JAVA_HOME%\bin;%FLUTTER_ROOT%\bin;%ANDROID_SDK_ROOT%\cmdline-tools\bin;%PATH%"
  
  echo =================================================================
  echo  环境已就绪!
  echo.
  echo  - Flutter 路径:     %FLUTTER_ROOT%
  echo  - Android SDK 路径: %ANDROID_SDK_ROOT%
  echo  - Java JDK 路径:    %JAVA_HOME%
  echo.
  echo  (如果设置了代理, 当前终端的所有网络请求将通过代理)
  echo.
  echo  正在启动 VS Code...
  echo =================================================================
  
  rem 启动 VS Code，它会继承当前这个终端的环境变量；使用 start异步启动 不阻塞当前终端 /B在当前窗口中启动程序，而不是创建一个新窗口
  start /B "" code .
  ```
## 其他配置
- flutter doctor 找不到 Android SDK
  ```
  flutter config --android-sdk "E:\Codes\Android\dev\android_sdk"
  ```
- 安装 Android SDK 组件，flutter可能会要求特定版本的Android SDK和build-tools，Android SDK一般不管，build-tools再下载一次需要版本即可（版本共存）
  ```
  sdkmanager --sdk_root="E:\Codes\Android\dev\android_sdk" --install "platform-tools" "platforms;android-34" "build-tools;34.0.0"
  ```
- 同意安卓许可证
  ```
  flutter doctor --android-licenses
  ```
- 最终检查
  ```
  flutter doctor
  ```
- Java Gradle 网络代理，%USERPROFILE%/.gradle下新建gradle.properties文件，内容如下：
  ```
  # ===================================================
  # == Gradle 全局代理配置
  # ===================================================
  
  # --- HTTP 代理 ---
  systemProp.http.proxyHost=127.0.0.1
  systemProp.http.proxyPort=10808
  systemProp.http.nonProxyHosts=*.nonproxyrepos.com|localhost
  
  # --- HTTPS 代理 ---
  systemProp.https.proxyHost=127.0.0.1
  systemProp.https.proxyPort=10809
  systemProp.https.nonProxyHosts=*.nonproxyrepos.com|localhost
  ```
# 常用命令
## gradle
- Gradle 锁文件问题
  ```
  taskkill /F /IM java.exe /T
  ```
- NDK缺少 source.properties，删除ndk目录下所有文件
  ```
  E:\Codes\Android\dev\android_sdk\ndk
  ```
## flutter
- 构建
  ```
  flutter run
  ```
- 清理构建
  ```
  flutter clean
  ```
- 
