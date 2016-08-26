# pycbg
Python implementation to geolocate IP addresses using SoL constraints.

Uses ping measurements from RIPE Atlas as input.
  
Example:
```python
cbgObj=cbg()
inputConstraints=cbgObj.getInputConstraints(json.load(open('data/4514096.json','r')))
cityDict=cbgObj.getCities(inputConstraints)
```