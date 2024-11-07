@echo off
set ssid=zzuli-student

echo Attempting to connect to WiFi network: %ssid%
netsh wlan connect name="%ssid%"
if %errorlevel% neq 0 (
    echo An error occurred while connecting to the WiFi network.
) else (
    echo Connected to WiFi network successfully.
)

REM 第一次获取并打印 IP 地址
for /f "tokens=2 delims=:" %%i in ('ipconfig ^| findstr /i "IPv4"') do set IP=%%i
set IP=%IP:~1%
if defined IP (
    echo First IP address: %IP%
) else (
    echo Unable to retrieve IP address.
)

REM 等待 5 秒
timeout /t 5 /nobreak >nul

REM 第二次获取并打印 IP 地址
for /f "tokens=2 delims=:" %%i in ('ipconfig ^| findstr /i "IPv4"') do set IP=%%i
set IP=%IP:~1%
if defined IP (
    echo Second IP address: %IP%
) else (
    echo Unable to retrieve IP address.
)

REM 准备登录请求参数
set student_account=332211020770
set student_account_pwd=2829992550Abc..
set url=http://10.168.6.10:801/eportal/portal/login

@REM REM 使用 curl 发送请求
@REM curl -G "%url%" ^
@REM     --data-urlencode "callback=dr1003" ^
@REM     --data-urlencode "login_method=1" ^
@REM     --data-urlencode "user_account=,0,%student_account%@zzulis" ^
@REM     --data-urlencode "user_password=%student_account_pwd%" ^
@REM     --data-urlencode "wlan_user_ip=%IP%" ^
@REM     --data-urlencode "wlan_user_ipv6=" ^
@REM     --data-urlencode "wlan_user_mac=000000000000" ^
@REM     --data-urlencode "wlan_ac_ip=10.168.6.9" ^
@REM     --data-urlencode "wlan_ac_name=" ^
@REM     --data-urlencode "jsVersion=4.2.1" ^
@REM     --data-urlencode "terminal_type=1" ^
@REM     --data-urlencode "lang=zh-cn" ^
@REM     --data-urlencode "v=2014"

@REM if %errorlevel% neq 0 (
@REM     echo An error occurred while sending the login request.
@REM ) else (
@REM     echo Login request sent successfully.
@REM )

@REM pause

REM 使用 curl 发送 GET 请求
set "student_account=332211020770"
set "student_account_pwd=2829992550Abc.."
set "url=http://10.168.6.10:801/eportal/portal/login?callback=dr1003&login_method=1&user_account=%%2C0%%2C%student_account%@zzulis&user_password=%student_account_pwd%&wlan_user_ip=%IP%&wlan_user_ipv6=&wlan_user_mac=000000000000&wlan_ac_ip=10.168.6.9&wlan_ac_name=&jsVersion=4.2.1&terminal_type=1&lang=zh-cn&v=2014"

curl "%url%"

if %errorlevel% neq 0 (
    echo An error occurred while sending the login request.
) else (
    echo Login request sent successfully.
)

pause
