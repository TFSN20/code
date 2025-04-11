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
- 设置代理。
  ```
  regexp:.*\.docker\.com.*,
  regexp:.*\.docker\.io.*
  ```
- 设置镜像文件夹位置：docker desktop -> settings -> Resources -> Advanced -> Disk image location(brower) -> Apply & restart(等待一会)
- 进入FastGPT目录，下载镜像
  ```
  docker-compose -f docker-compose-pgvector.yml pull
  ```
  如果yml名为docker-compose.yml，可以直接
  ```
  docker-compose pull
  ```
## 使用
- 启动容器
  ```
  docker-compose -f docker-compose-pgvector.yml up -d
  ```
  如果yml名为docker-compose.yml，可以直接
  ```
  docker-compose up -d
  ```
- 访问fastgpt
  ```
  localhost:3000
  ```
- 最新版默认使用AIProxy进行模型api统一化，如果想要之前版本的one api则需要到yml文件夹中修改。这里我们使用AIProxy，集成到了FastGpt中。
  账号 -> 模型提供商 -> 模型渠道
  - 模型渠道
    - 定义常见模型渠道（如openai，gemini，豆包等）的统一key和api（api有默认值，可修改）。
    - 模型ID十分重要，这是请求中区分使用的哪一个模型的标识。
    - 渠道名仅用于区分，可添加内置的模型，也可自定义添加模型（模型ID保持一致），这样我们可以将需要的模型添加进来，或自定义模型（比如一些较新的）
  ```
  localhost:3001
  ```
## 常用api
- gemini，好啊好啊， 1,000,000 TPM + 1,500 RPD + 15 RPM，用这个分析小说简直不要太美。
  ```
  https://aistudio.google.com/apikey
  ```
- chateverything，免费，支持openai索引模型和3.5turbo，但input限制4096tokens
  ```
  https://github.com/chatanywhere/GPT_API_free
  ```
- 硅基流动，注册就送14元，邀请一人得14元。美。fastgpt对于硅基流动没有内置太多模型，更多的需要自己手动添加模型ID，或进行模型映射（名称可以点击模型直接复制的哪个，通常带有/）。
  ```
  https://cloud.siliconflow.cn/account/ak
  ```
- 同义千问，100万tokens
  ```
  https://bailian.console.aliyun.com/?tab=model#/api-key
  ```
- openai不建议使用，扣扣索索，还会过期。grok需要消费5美元，并开启数据共享，每月就能享受150美元免费额度。
## 常用命令
- 查看镜像
  ```
  docker images
  ```
- 拉取镜像
  ```
  docker-compose pull
  ```
- 启动容器
  ```
  docker-compose up -d
  ```
- 容器详细
  ```
  docker ps
  ```
## 常用功能
- 容器管理，镜像查看，容器详细都可以在docker desptop进行。
