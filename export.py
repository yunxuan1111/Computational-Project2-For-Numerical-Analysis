import comp2
from comp2 import run
from pathlib import Path  

df,df_2 = run()


filepath = Path('D:\Comp 2\MethodIntegrand.csv')  
filepath.parent.mkdir(parents=True, exist_ok=True)  
filepath2 = Path('D:\Comp 2\IntegrandMethod.csv')  
filepath2.parent.mkdir(parents=True, exist_ok=True)  
df.to_csv(filepath)
df_2.to_csv(filepath2)

