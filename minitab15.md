# 箱线图
## 通用命令
- 文本使用双引号
  ```
  ""1""
  ```
- 结束命令用英文句号
## 两个分类列的箱线图
- Ctrl + L调出命令行编辑器
- 一个Boxplot一个图
  ```
  Boxplot '揉平后正极耳外漏' * '磁悬浮动子';
    Include;
      Where "'揉平工位' = 1";
    Title "揉平工位 = 1 的箱线图";
    Reference 2 0.2 1.6;
      Color 2;
      Type 2;
    Reference 2 0.7 0.9 1.1;
      Color 2;
      Type 2.
  Boxplot '揉平后正极耳外漏' * '磁悬浮动子';
    Include;
      Where "'揉平工位' = 2";
    Title "揉平工位 = 1 的箱线图";
    Reference 2 0.2 0.7 0.9 1.1 1.6;
      Color 2;
      Type 2.
  Boxplot '揉平后正极耳外漏' * '磁悬浮动子';
    Include;
      Where "'揉平工位' = 3";
    Title "揉平工位 = 1 的箱线图";
    Reference 2 0.2 0.7 0.9 1.1 1.6;
      Color 2;
      Type 2.
  Boxplot '揉平后正极耳外漏' * '磁悬浮动子';
    Include;
      Where "'揉平工位' = 4";
    Title "揉平工位 = 1 的箱线图";
    Reference 2 0.2 0.7 0.9 1.1 1.6;
      Color 2;
      Type 2.
  Boxplot '揉平后正极耳外漏' * '磁悬浮动子';
    Include;
      Where "'揉平工位' = 5";
    Title "揉平工位 = 1 的箱线图";
    Reference 2 0.2 0.7 0.9 1.1 1.6;
      Color 2;
      Type 2.
  Boxplot '揉平后正极耳外漏' * '磁悬浮动子';
    Include;
      Where "'揉平工位' = 6";
    Title "揉平工位 = 1 的箱线图";
    Reference 2 0.2 0.7 0.9 1.1 1.6;
      Color 2;
      Type 2.
  Boxplot '揉平后正极耳外漏' * '磁悬浮动子';
    Include;
      Where "'揉平工位' = 7";
    Title "揉平工位 = 1 的箱线图";
    Reference 2 0.2 0.7 0.9 1.1 1.6;
      Color 2;
      Type 2.
  Boxplot '揉平后正极耳外漏' * '磁悬浮动子';
    Include;
      Where "'揉平工位' = 8";
    Title "揉平工位 = 1 的箱线图";
    Reference 2 0.2 0.7 0.9 1.1 1.6;
      Color 2;
      Type 2.
  ```
