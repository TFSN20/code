- 下载并安装Android Studio
- 安装时有Sdk选项，SDK全称为Software Development Kit,意即软件开发工具包,它可将App源码编译为可执行的App应用。
- 新版AS新建Empty Activity项目时只能选择Kotlin语言开发。
- 如果之前安装过，则十分建议删除用户名目录下的.android和.gradle文件夹，AS的安装目录，Sdk所在的目录，以及项目目录（很重要）。

## gradle问题
安装会在用户名目录下创建.android和.gradle文件夹，其中.gradle\wrapper\dists下存放各种版本的gradle，对应gradle-x.x-bin目录，
每个目录下的一连串数字字母符号下的gradle-x.x文件夹才是bin存放的地方，其中gradle-x.x同级可能存在存在（需要Sync刷新）gradle-x.x-bin.zip.lck和gradle-x.x-bin.zip.ok，gradle-x.x-bin.zip.part文件夹。
AS中gradle版本下载很慢，可以直接把外部下载好的zip放到一连串数字字母符号文件夹下，并解压（也就是说其实zip和解压文件都要有，但是好像这个方法没用，还是用本地代理完成全部下载吧）。
- 参考：https://blog.csdn.net/u011046452/article/details/107529346

下载完gradle还会下载maven，可以使用本地代理，设置127.0.0.1，v2ray中端口（分为走http和socks端口不一样）是10809，输入网址（https和http不一样）检查是否走了v2ray(然后还要重启AS，会再探出一个代理窗口，用于http和https)。
- 参考：https://blog.csdn.net/qq_65732918/article/details/134903527

## 调试
USB真机调试：设置——搜索Sdk——SDK Tools——勾选并点击下载Google USB Driver——Android版本要和手机Android版本一致——手机打开开发者调试和安装。
参考：https://www.cnblogs.com/rainbow70626/p/14364552.html，https://cloud.tencent.com/developer/article/1743279
