- https://mp.weixin.qq.com/
- 注册小程序：未注册过任何公众号或小程序的邮箱
- 开发与服务-开发管理-AppID

## 开发
- 选用vscode，安装wxml插件
- 规定div为view
### 云函数
- 目录结构:  
  cloudfunctions/  
    └── quickstartFunctions/  
        ├── getOpenId/  
        │   └── index.js      // getOpenId 逻辑  
        ├── changeVote.js     // changeVote 逻辑  
        ├── index.js          // 云函数入口文件  
        └── package.json  
- index.js
  ```
  const getOpenId = require('./getOpenId/index');
  const changeVote = require('./changeVote.js');
  // 云函数入口函数
  exports.main = async (event, context) => {
    switch (event.type) {
      case 'getOpenId':
        return await getOpenId.main(event, context);
        case 'changeVote':
          return await changeVote(event, context);
    }
  };
  ```
- getOpenId：index.js
  ```
  const cloud = require('wx-server-sdk');
  cloud.init({
    env: cloud.DYNAMIC_CURRENT_ENV
  });
  // 获取openId云函数入口函数
  exports.main = async (event, context) => {
    // 获取基础信息
    const wxContext = cloud.getWXContext();
    return {
      openid: wxContext.OPENID,
      appid: wxContext.APPID,
      unionid: wxContext.UNIONID,
    };
  };
  ```
- changeVote.js：
  ```
  const cloud = require('wx-server-sdk');
  cloud.init({
    env: 'tfsn20-qqt-1gt71yjv8aqqf8d282'
  });
  
  const db = cloud.database();
  
  module.exports = async (event, context) => {
    const productId = event.productId; // 从前端传来的产品ID
    const priceIndex = event.priceIndex;
    return {}
  }
  ```
- 调用云函数之前先进行上传：
  ```
  wx.cloud.callFunction({
    name: 'quickstartFunctions',
    data: {
      type: 'changeVote',  // 指定要调用的子函数
      productId: '25essssssxxxxx9f782f8',
      priceIndex:this.data.selectedOptionIndex
    },
    success: res => {
      if (res.result.success) {
        wx.showToast({
          title: res.result.message,
        });
      } else {
        wx.showToast({
          title: res.result.message,
          icon: 'none'
        });
      }
    },
    fail: err => {
      wx.showToast({
        title: '购买失败，请重试',
        icon: 'none'
      });
      console.error('调用云函数失败', err);
    }
  });
  ```
- 目录与文件的解析规则:  
  - 当 require 的路径指向一个目录（如 ./getOpenId），Node.js 会尝试读取该目录下的 index.js 文件。
    这是因为 index.js 是一个目录默认入口文件，所以你可以写 require('./getOpenId') 或require('./getOpenId/index')，而不需要添加 .js 后缀。  
  - 当 require 的路径指向具体的文件（如 ./changeVote.js），如果省略 .js，Node.js 会尝试自动补全文件扩展名，但它会先尝试 .json、.node 等其他扩展名，再找到 .js 文件，因此并不总是能成功解析。
    所以直接写 require('./changeVote.js') 更加保险。
  
  

  

