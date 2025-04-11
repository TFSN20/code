## 下载
- 下载Ubuntu。只能安装在C盘，可使用mklink更改位置。
- 管理员cmd执行以下命令，将wsl更改为2, Ubuntu-22.04为下载的Ubuntu。
  ```
  wsl --list --verbose
  ```
  ```
  wsl --set-version Ubuntu-22.04 2
  ```
- 下载docker桌面版。下载完需要重启电脑。
- 新建一个目录名称任意，这里为FastGPT，放入config.json和docker-compose-pgvector.yml（docker-compose-xxx.yml都行，个人娱乐一般用pgvector）文件。
- docker设置里的Docker Engine json配置
  ```
  {
  "builder": {
    "gc": {
      "defaultKeepStorage": "20GB",
      "enabled": true
    }
  },
  "experimental": false,
  "registry-mirrors": [
    "https://registry.ggnb.top"
  ]
  }
  ```
- 消注释阿里云镜像下载源image行，注释掉原来的下载源。
- 进入FastGPT目录，下载镜像
  ```
  docker-compose -f docker-compose-pgvector.yml pull
  ```
  如果yml名为docker-compose.yml，可以直接
  ```
  docker-compose pull
  ```
