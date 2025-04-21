## MCP
- Model Context Protocol (模型上下文协议)
- AI：LLMs（大语言模型）
- MCP server：AI与外部工具的中间层。
  - 每个MCP服务都是一个工具，每个MCP服务通常是一段nodejs代码或一个python程序。
  - AI通过操作系统的stdio（标准输入通道）调用某个MCP server，消息格式可以是json，MCP server接收到消息后执行自己的代码功能（nodejs代码或python程序）或使用API请求访问外部工具。
- 过程详解：
  - 用户输入问题，然后支持MCP协议的客户端（MCP client）和AI进行对接，主流对接有两种，系统提示词和function call。
    - claude：这个支持MCP协议的客户端发送的信息包括system系统提示词的关于MCP的使用，以及有哪些MCP server及其参数。
    - 5ire：这个支持MCP协议的客户端则是使用tool字段，需要AI支持function call功能。
  - AI返回固定的调用MCP的格式消息，MCP client根据这个消息直接调用对应MCP server，并将MCP server固定的调用MCP的格式消息加上用户的问题再次发送到AI，最终用户看到的是此次AI的结果。
![image](https://github.com/user-attachments/assets/2b42bfab-43d3-4ec7-8d7e-33c7efad7bd9)

