# Coding Specifications


## 1. 注释的规范(类似JavaDoc)
计划注释的使用格式
```
/**
 * @brief [what this function does]
 * @note [some instruction and explanation]
 * @param arg1
 * @param arg2
 * ...
 * @retVal/@return ...
 */
int func(int arg1, char *arg2, MyStruct *arg3, ...){
  /*
   * Things that cannnot be explained by codes, expecially
   * some static or local variables, and the reasons of why coding
   * the function like this. 
   */
  [procedure 1]
  // partrition the code for better reading and designing
  [procedure 2]
  return value;
}
```
## 2. 代码的规范(Google风格)
```
RetType funcTest(ArgType1 arg1, ArgType2 arg2, ...) {
  LocalVarType1 var1;
  LocalVarType2 varShortName;
  LocalVarType3 varFunc_LongName;
  // ...
  
  for(int i = 0; i < argX; i++) {
    // ...
    if( ... ) {
      // ...
    } else {
      // ...
    }
    ...
  }
  
  // ...
  funcAnother(varR, varQ, ...);
}
```

## 3. 编译运行的规范

编译指令：gcc *.c -o main -lhts

程序使用：./main \<command\> \[arguments\]

Usage可以直接通过命令行"./main"显示