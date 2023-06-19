mod_biome <- function(Z) 
{
  npp=Z[3]
  pft=Z[15:27]
  optpft=Z[4]
  gdd=Z[10]
  mtco=Z[6]
  tcm=mtco-15
  alpha=Z[8]
  biome=Z[2]
  
  
  dpft=c(8:10,14)
  
  # steppes / alpha
  if (alpha<=33 & gdd>300) {
    if (gdd>6000 | tcm>0) {biome=19} else {biome=20} 
  }
  # desert / alpha
  if (alpha<=18 & gdd>300) {
    if (gdd>6000 | tcm>0) {biome=21} else {biome=18} 
  }
  
  # steppes ou savanne
  if (any(optpft==dpft) & alpha<=45) {
    if (gdd>6000 | tcm>0) {biome=19} else {biome=20} 
    if ((gdd>6000 | tcm>0) & optpft==14) {biome=12} 
  }
  
  # desert
  if (any(optpft==dpft) & alpha<=28) {
    if (gdd>6000 | tcm>0) {biome=21} else {biome=18} 
  }
  
  if (gdd<10) {biome=27}
  
  return (biome)
}
