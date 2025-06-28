- 不同文件后缀的编码, ctrl + shift + p 输入 settings.json ，打开用户设置
  ```
    "files.associations": {
        "*.ps1": "powershell"
    },
    // 第二部分：为 "powershell" 语言设置特定的编码
    // 作用：告诉 VS Code，对于所有被识别为 PowerShell 语言的文件，
    // 在保存时都使用 utf8bom 编码，这会覆盖全局的默认编码(utf8)。
    "[powershell]": {
        "files.encoding": "utf8bom"
    }
  ```
  
