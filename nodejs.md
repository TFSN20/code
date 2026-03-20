## 安装
- 可以安装到D盘
- 查看版本
  ```
  node -v
  npm -v
  ```
- 全局安装位置修改
  - 在nodejs安装目录下新建两个目录
    ```
    node_global
    ```
    ```
    node_cache
    ```
  - cmd执行
    ```
    npm config set prefix "node_global路径"
    ```
    ```
    npm config set cache "node_cache路径"
    ```
  - 修改全局模块prefix位置的用户环境变量
    ```
    C:\Users\Administrator\AppData\Roaming\npm 改为 D:\ProgramFiles\Codes\nodejs\node_global
    ```
  - 新增NODE_PATH系统变量
    ```
    D:\ProgramFiles\Codes\nodejs\node_global\node_modules
    ```
  - 修改镜像源（默认https://registry.npmjs.org/）
    ```
    npm config set registry https://registry.npmmirror.com/
    ```
  - 文件夹权限
    ```
    赋予Authenticated Users用户完全控制权限对于nojs安装目录
    ```
- 常用功能
  - 查看模块安装位置
    ```
    npm root
    ```
  - 查看模块列表
    ```
    npm list
    ```
  - 查看模块顶层依赖
    ```
    npm list --depth=0
    ```
  - 查看模块缓存位置
    ```
    npm config get cache
    ```
  - 查看全局模块prefix
    ```
    npm config get prefix
    ```
  - 查看下载源
    ```
    npm config get registry
    ```
  - 查看配置文件位置
    ```
    npm config get userconfig
    ```
  - 初始化项目
    ```
    npm init
    ```
    ```
    npm init -y
    ```
  ## fnm安装和使用
  - fnm.exe github下载 并将下载目录添加到用户Path变量
  - fnm 环境变量信息
    ```
    fnm env
    ```
  - node安装位置修改：新建指定文件夹，将路径添加到用户变量，名为`FNM_DIR`
  - node下载源修改：用户变量新建`FNM_NODE_DIST_MIRROR	`，值为`https://npmmirror.com/mirrors/node`
  - 查看可安装版本
    ```
    fnm ls-remote
    ```
  - 安装最新长期稳定版本
    ```
    fnm install --lts
    ```
  - 安装指定版本最新稳定版（安装 22.x 的最新稳定版）
    ```
    fnm install 22
    ```
  - 加载node文件夹（仅当前powershell窗口）
    ```
    fnm env --use-on-cd | Out-String | Invoke-Expression
    ```
  - 设置全局默认版本
    ```
    fnm default 22
    ```
  - 查看已安装版本
    ```
    fnm ls
    ```
    
