----
软链命令
ln -s /SGRNJ06/randd/USER/wangjingshen/soft/ossutil/ossutil-v1.7.19-linux-amd64/ossutil64


----
step1.ossutil配置

ossutil64 config
#Please enter endpoint:         输入  oss-cn-hangzhou.aliyuncs.com
#Please enter stsToken:         直接回车
#Please enter accessKeyID:      输入ID
#Please enter accessKeySecret:  输入Secret


step2.下载
ossutil64 cp -r oss://examplebucket/destfolder/ localfolder/

例如: ossutil64 cp -r oss://singleronbio-data-release/PK22011902_GueJ/ test/

