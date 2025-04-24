- 需要下载额外才能抓取
- 过滤示例
  ```
  (string(ipv6.addr) matches "2409:8c00:.*" or string(ipv6.addr) matches "2409:8c20:.*" or string(ipv6.addr) matches "2409:8c44:.*" or string(ip.addr) matches "36.156..*" or string(ip.addr) matches "101.42..*" or string(ip.addr) matches "39.156..*" or string(ip.addr) matches "175.102..*")
  ```
  ```
  http and http.request.uri contains "playUrl" 
  ```
