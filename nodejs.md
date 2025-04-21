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
    
