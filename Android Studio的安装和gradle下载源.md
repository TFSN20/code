## 安装
- 如果之前安装过，则十分建议删除用户名目录下的.android和.gradle文件夹，AS的安装目录，Sdk所在的目录，以及项目目录（很重要）。
- 下载并安装Android Studio，安装会在用户名目录下创建.android和.gradle文件夹。
- 安装时有Sdk选项，SDK全称为Software Development Kit,意即软件开发工具包,它可将App源码编译为可执行的App应用。
## gradle问题
- 安装会在用户名目录下创建.android和.gradle文件夹，其中.gradle\wrapper\dists下存放各种版本的gradle，对应gradle-x.x-bin目录，
每个目录下的一连串数字字母符号下的gradle-x.x文件夹才是bin存放的地方，其中gradle-x.x同级可能存在存在（需要Sync刷新）gradle-x.x-bin.zip.lck和gradle-x.x-bin.zip.ok，gradle-x.x-bin.zip.part文件夹。
- AS中gradle版本下载很慢，可以直接把外部下载好的zip放到一连串数字字母符号文件夹下，并解压（也就是说其实zip和解压文件都要有，但是好像这个方法没用，还是用本地代理完成全部下载吧）。
- 参考：https://blog.csdn.net/u011046452/article/details/107529346
- 下载完gradle还会下载maven，可以使用本地代理，设置127.0.0.1，v2ray中端口（分为走http和socks端口不一样）是10809，输入网址（https和http不一样）检查是否走了v2ray(然后还要重启AS，会再探出一个代理窗口，用于http和https)。
- 参考：https://blog.csdn.net/qq_65732918/article/details/134903527
## 调试
- USB真机调试：设置——搜索Sdk——SDK Tools——勾选并点击下载Google USB Driver——Android版本要和手机Android版本一致——手机打开开发者调试和安装。
参考：https://www.cnblogs.com/rainbow70626/p/14364552.html，https://cloud.tencent.com/developer/article/1743279
## 新建项目
- 新版AS新建Empty Activity项目时只能选择Kotlin语言开发。
- 新建java文件：
  以Project结构查看项目结构，在app/src/main/java/包名下的AB.java文件（AS中不显示java后缀）和app/src/main/res/layout下的b_a.xml是一对，表示逻辑和页面，新建java时（就是Activity文件）对应的xml会自动创建。
- 每个xml相互独立，因此可以有相同的ID的Button。
- 如果你想要在一个Activity中操作另一个Activity的布局元素，你通常需要采取间接的方式，例如：
通过Intent传递数据：在AActivity中启动BActivity时，可以通过Intent传递信息给BActivity，告诉它如何更新自己的界面。
使用全局变量或单例：通过一个共享的数据存储（如应用的全局状态或单例模式类），在不同的Activity之间同步状态或数据。
更新数据源：如果两个Activity依赖同一个数据源（例如数据库），那么一个Activity更新了数据源后，另一个Activity可以根据这些变更来更新其界面。
- 每个页面都必须有个盒子
- xml添加元素
  ```
      <Button
        android:id="@+id/uploadBtn"
        android:layout_width="wrap_content"
        android:layout_height="wrap_content"
        app:layout_constraintBottom_toBottomOf="parent"
        app:layout_constraintTop_toTopOf="parent"
        app:layout_constraintStart_toStartOf="parent"
        app:layout_constraintEnd_toEndOf="parent"
        android:text="上传"
        />
  ```
- Activity.java中查找元素
  ```
  Button recommendBtn = findViewById(R.id.recommendBtn);
  ```
- 切换页面（Intent）
  ```
  Intent uploadIntent = new Intent(MainActivity.this, UploadActivity.class);
  startActivity(uploadIntent);
  ```

