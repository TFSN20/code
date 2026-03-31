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
- 8个图组合为两个组图
  ```
  Layout;
    Title "412B 1#~4# 揉平工位".
  Boxplot '揉平后正极耳外漏' * '磁悬浮动子';
    Include;
      Where "'揉平工位' = 1";
    Title "揉平工位 = 1 的箱线图";
    Reference 2 0.2 1.6;
      Color 2;
      Type 2;
    Reference 2 0.7 0.9 1.1;
      Color 2;
      Type 2;
    Figure 0.0 0.5 0.5 1.0.
  Boxplot '揉平后正极耳外漏' * '磁悬浮动子';
    Include;
      Where "'揉平工位' = 2";
    Title "揉平工位 = 2 的箱线图";
    Reference 2 0.2 0.7 0.9 1.1 1.6;
      Color 2;
      Type 2;
    Figure 0.5 1.0 0.5 1.0.
  Boxplot '揉平后正极耳外漏' * '磁悬浮动子';
    Include;
      Where "'揉平工位' = 3";
    Title "揉平工位 = 3 的箱线图";
    Reference 2 0.2 0.7 0.9 1.1 1.6;
      Color 2;
      Type 2;
    Figure 0.0 0.5 0.0 0.5.
  Boxplot '揉平后正极耳外漏' * '磁悬浮动子';
    Include;
      Where "'揉平工位' = 4";
    Title "揉平工位 = 4 的箱线图";
    Reference 2 0.2 0.7 0.9 1.1 1.6;
      Color 2;
      Type 2;
    Figure 0.5 1.0 0.0 0.5.
  EndLayout.
  
  Layout;
    Title "412B 5#~8# 揉平工位".
  Boxplot '揉平后正极耳外漏' * '磁悬浮动子';
    Include;
      Where "'揉平工位' = 5";
    Title "揉平工位 = 5 的箱线图";
    Reference 2 0.2 0.7 0.9 1.1 1.6;
      Color 2;
      Type 2;
    Figure 0.0 0.5 0.5 1.0.
  Boxplot '揉平后正极耳外漏' * '磁悬浮动子';
    Include;
      Where "'揉平工位' = 6";
    Title "揉平工位 = 6 的箱线图";
    Reference 2 0.2 0.7 0.9 1.1 1.6;
      Color 2;
      Type 2;
    Figure 0.5 1.0 0.5 1.0.
  Boxplot '揉平后正极耳外漏' * '磁悬浮动子';
    Include;
      Where "'揉平工位' = 7";
    Title "揉平工位 = 7 的箱线图";
    Reference 2 0.2 0.7 0.9 1.1 1.6;
      Color 2;
      Type 2;
    Figure 0.0 0.5 0.0 0.5.
  Boxplot '揉平后正极耳外漏' * '磁悬浮动子';
    Include;
      Where "'揉平工位' = 8";
    Title "揉平工位 = 8 的箱线图";
    Reference 2 0.2 0.7 0.9 1.1 1.6;
      Color 2;
      Type 2;
    Figure 0.5 1.0 0.0 0.5.
  EndLayout.
  ```
# DOE
## 全因子设计
- 常用1-4个因子，具体操作：
  ```
  2水平因子（默认生成元） -> 因子数 -> 设计 ->
                                            全因子 -> 每个区域的中心点数（中心点实验次数） -> 角点的仿行数（一般为1 即只实验一次 无重复 ） -> 区组数（一般为1，即同一天 同一设备 同一环境之意）
                                         -> 因子（设置高低值）
  ```
- 分析
  ```
  分析因子设计 -> 项 -> 所选项里添加或删除知道满足需求 ->
                 图形 -> 勾选正态和Pareto、四合一 -> 
                 结果 -> 一般选择系数和方差分析表(C)
  ```
    
