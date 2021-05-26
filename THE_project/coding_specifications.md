# Coding Specifications


## 1. Comments (JavaDoc style)

```
/**
 * @brief [what this function does]
 * @note [some instruction and explanation. Mainly personal blahs]
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
  // ---------------------- procedure 1 -----------------------------
  [procedure 1]
  // ---------------------- procedure 2 -----------------------------
  [procedure 2]
  return value;
}
```
## 2. Codes (Google style)
```
RetType funcTest(ArgType1 arg1, ArgType2 arg2, ...) {
  LocalVarType1 var1;
  LocalVarType2 varShortName;
  LocalVarType3 varFunc_LongName;
  // ...
  
  for(int i = 0; i < argX; i++) {
    LocalVarType4 tmp_Var;
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

## 3. Naming Rules
```
1. type_operation
2. typeOperation
```