## 版本
- Auto.js_4.1.1 Alpha2：支持importClass，不支持模板字符串
- Autojsx：最新版不支持importClass，支持模板字符串。
## 常见问题
- http模块请求时最好不携带accept-encoding请求头，否则响应体大概率乱码。
- 关于悬浮窗和UI上的界面操作（按钮点击），如果有耗时操作会堵塞，使其不能点击，使用threads.start(()=>{})解决。
- jsBridge.callHandler里的网络请求必须threads.start。
- 悬浮窗按钮点击后触发的函数若有堵塞，则必须先 yield myAwait(sleep_(10));
- 关于UI的操作必须在ui.run里进行，且上层函数必须是ui.run
- 响应体"Accept-Encoding", "deflate"。string乱码
  ```
  importClass(java.io.ByteArrayOutputStream);
  importClass(java.io.InputStreamReader);
  importClass(java.io.BufferedReader);
  importClass(java.lang.StringBuffer);
  importClass(java.net.HttpURLConnection);
  importClass(java.net.URL);
  importClass(java.util.zip.Inflater);
  importClass(java.util.zip.InflaterInputStream);
  
  function decompressDeflate(inputStream) {
      let inflater = new Inflater(true); // true 表示跳过 zlib 头
      let buffer = java.lang.reflect.Array.newInstance(java.lang.Byte.TYPE, 1024);
      let byteArrayOutputStream = new ByteArrayOutputStream();
  
      let len;
      while ((len = inputStream.read(buffer)) !== -1) {
          byteArrayOutputStream.write(buffer, 0, len);
      }
  
      let decompressedBytes = byteArrayOutputStream.toByteArray();
      byteArrayOutputStream.close();
  
      inflater.setInput(decompressedBytes);
      let outputStream = new ByteArrayOutputStream();
      let outBuffer = java.lang.reflect.Array.newInstance(java.lang.Byte.TYPE, 1024);
  
      try {
          while (!inflater.finished()) {
              let count = inflater.inflate(outBuffer);
              outputStream.write(outBuffer, 0, count);
          }
      } finally {
          inflater.end();
      }
  
      return outputStream.toString("UTF-8");
  }
  
  function fetchBilibiliDM(urlPath) {
      let conn = (new URL(urlPath)).openConnection();
      conn.setConnectTimeout(5000);
      conn.setRequestMethod("GET");
      conn.setRequestProperty("Accept-Encoding", "deflate");
  
      conn.connect();
  
      if (conn.getResponseCode() === 200) {
          let inputStream = conn.getInputStream();
          let responseString = decompressDeflate(inputStream);
  
          console.log("原始响应体字符串:");
          console.log(responseString);
  
          return responseString;
      } else {
          console.error("HTTP 请求失败，状态码: " + conn.getResponseCode());
          return null;
      }
  }
  
  // 示例调用
  let cid = "12345678"; // 替换为实际的 cid
  let url = `https://api.bilibili.com/x/v1/dm/list.so?oid=${cid}`;
  let response = fetchBilibiliDM(url);
  if (response) {
      console.log("成功获取响应体字符串:");
      console.log(response);
  } else {
      console.error("请求失败或无效响应。");
  }
  ```
  - java.util.zip.ZipException: incorrect header check 通常是由于解压方式与数据的实际压缩方式不匹配导致的。在你提供的示例中，我们假设响应的 Content-Encoding 是 deflate，但如果服务器返回的压缩格式并不完全符合标准（例如，部分服务会返回原始的 zlib 数据头），可能会导致这个错误。若返回的数据是原始 zlib 数据头而非标准的 deflate，可以尝试直接使用 Inflater 而不是 InflaterInputStream 进行解压：
    ```
    importClass(java.io.ByteArrayOutputStream);
    importClass(java.io.InputStreamReader);
    importClass(java.io.BufferedReader);
    importClass(java.lang.StringBuffer);
    importClass(java.net.HttpURLConnection);
    importClass(java.net.URL);
    importClass(java.util.zip.Inflater);
    importClass(java.util.zip.InflaterInputStream);
    
    function decompressDeflate(inputStream) {
    let inflater = new Inflater(true); // true 表示跳过 zlib 头
    let buffer = java.lang.reflect.Array.newInstance(java.lang.Byte.TYPE, 1024);
    let byteArrayOutputStream = new ByteArrayOutputStream();

    let len;
    while ((len = inputStream.read(buffer)) !== -1) {
        byteArrayOutputStream.write(buffer, 0, len);
    }

    let decompressedBytes = byteArrayOutputStream.toByteArray();
    byteArrayOutputStream.close();

    inflater.setInput(decompressedBytes);
    let outputStream = new ByteArrayOutputStream();
    let outBuffer = java.lang.reflect.Array.newInstance(java.lang.Byte.TYPE, 1024);

    try {
        while (!inflater.finished()) {
            let count = inflater.inflate(outBuffer);
            outputStream.write(outBuffer, 0, count);
        }
    } finally {
        inflater.end();
    }

    return outputStream.toString("UTF-8");
    }
    
    function fetchBilibiliDM(urlPath) {
    let conn = (new URL(urlPath)).openConnection();
    conn.setConnectTimeout(5000);
    conn.setRequestMethod("GET");
    conn.setRequestProperty("Accept-Encoding", "deflate");

    conn.connect();

    if (conn.getResponseCode() === 200) {
        let inputStream = conn.getInputStream();
        let responseString = decompressDeflate(inputStream);
        return responseString;
    } else {
        console.error("HTTP 请求失败，状态码: " + conn.getResponseCode());
        return null;
    }
    }
    
    // 示例调用
    let cid = "28059239264"; // 替换为实际的 cid
    let response = fetchBilibiliDM(`https://api.bilibili.com/x/v1/dm/list.so?oid=${cid}`);
    ```
