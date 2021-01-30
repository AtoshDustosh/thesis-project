# Stage 2 Objectives

## 一些规范

#### 1. 注释的规范
计划注释的使用格式
```
/**
 * @brief [what this function does]
 * @note [some instruction and explanation]
 * @param arg1
 * @param arg2
 * ...
 * @return ...
 */
int func(int arg1, char *arg2, MyStruct *arg3, ...){
  /*
   * Things that cannnot be explained by codes, expecially
   * some static or local variables.
   */
  [procedure 1]
  // partrition the code for better reading and designing
  [procedure 2]
  return value;
}
```
#### 2. 代码的规范
```
RetType funcTest(ArgType1 arg1, ArgType2 arg2, ...)
{
  LocalVarType1 var1;
  LocalVarType2 varShortName;
  LocalVarType3 varFunc_LongName;
  // ...
  
  for(int i = 0; i < argX; i++)
  {
    // ...
    if( ... )
    {
      // ...
    }
    else
    {
      // ...
    }
    ...
  }
  
  // ...
  funcAnother(varR, varQ, ...);
}
```

#### 3. 编译运行的规范

编译指令：gcc *.c -o main -lhts

程序使用：./main \<command\> \[arguments\]

for example, ./main -i data/example.vcf -g --firstlines 100，
表示设定输入文件为data/example.vcf，控制台输出该文件的前100行内容

### 4. 项目组织的规范

参见程序本身所在文件夹的文件组织