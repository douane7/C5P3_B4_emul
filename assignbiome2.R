assignbiome <- function(Z) 
  # aggregated pft  (8)
{
#PFT
#    1 = trt Tropical Trees 
#    2 = tbe Temperate Broadleaved Evergreen Trees 
#    3 = tet Temperate Trees
#    4 = bot Boreal Trees     
#    5 = teg C3/C4 temperate grass plant type
#    6 = trg4 C4 tropical grass plant type
#    7= wds C3/C4 woody desert plant type
#    8= tug Tundra grass
#Biomes
#1    1  Tropical evergreen forest
#2    2  Tropical semi-deciduous forest
#3    3  Tropical deciduous forest/woodland
#4    4  Temperate deciduous forest
#5    5  Temperate conifer forest
#6    6  Warm mixed forest
#7    7  Cool mixed forest
#7    8  Cool conifer forest
#11   9  Cold mixed forest
#11   10 Evegreen taiga/montane forest
#11   11 Deciduous taiga/montane forest
#12   12 Tropical savanna
#12   13 Tropical xerophytic shrubland
#14   14 Temperate xerophytic shrubland
#14   15 Temperate sclerophyll woodland
#14   16 Temperate broadleaved savanna
#11   17 Open conifer woodland
#11   18 Boreal parkland
#19   19 Tropical grassland
#20   20 Temperate grassland
#21   21 Desert
#23   22 Steppe tundra
#23   23 Shrub tundra
#23   24 Dwarf shrub tundra
#23   25 Prostrate shrub tundra
#23   26 Cushion forb lichen moss tundra
#27   27 Barren
#27   28 Land ice
#-------------------------------------------------------

npp=Z['NPPtot']
pft=Z[2:9]
optpft=which.max(pft)
gdd=Z['GDD5']
mtco=Z['MTCO']
tcm=mtco-15
biome=27

# barren & tundra
if (optpft==8) {
  if (mtco>15 & npp>1000) {biome=12}
  if (mtco>15 & npp<=1000) {biome=19}
  if (mtco<=15) {biome=20}
  if (gdd<200) {biome=23} 
  if (npp<10) {biome=27}
}

# desert
if (optpft==6) {biome=19}
if (optpft==7 | optpft==6) {
  if (mtco>15) {biome=13} 
  if (mtco<=15) {biome=14}
  if (npp<=100) {biome=21}
}

# temperate grass
if (optpft==5) {
  if (gdd<150) {biome=22} else {biome=21}
  if (mtco>15) {biome=19}
}

# bot
if (optpft==4) {
  biome=9
  if (gdd>900 & tcm> -19) {biome=5}
  if (gdd>900 & tcm<= -19) {biome=10}
  if (mtco>15) {biome=1}
} 

# tbe
if (optpft==2) { 
  if(gdd>3000) {biome=6} else {biome=7}
}

# tet
if (optpft==3) {
  if (tcm<=-15) {biome=9} else {biome=7}
  if (gdd>=3000 & tcm>3) {biome=6}else {biome=4}
}


# tropical forest
if (optpft==1) {biome=1}

# savanna
if (any(optpft==5:8) & mtco>=15 & npp>1000) {biome=12}
if (any(optpft==5:8) & mtco<15 & npp>1000) {biome=14}
if (any(optpft==5:8) & gdd<900 & npp>1000) {biome=18}


return (biome)
}



