# load the dll
path_to_dll <- "C:\\Users\\jb240893\\Downloads\\vol-lib-dll-20250701\\VolLibDll20250701\\vollib64" # change this for your machine as needed
dyn.load(file.path(path_to_dll,"vollib.dll"))

# get version number of Fortran (just fyi)
.Fortran("vernum_r",vernum=integer(1))

# 1 retrieve default equation number for some species and area

# set some basic input information, using specific variable types for 
#  compatiblity with fortran
region <- as.integer(1)      # northern usfs region
forest <- as.character(16)   # lolo national forest
district <- as.character("01") #
species <- as.integer(202)   # Douglas-fir

# define the variables that are RETURNED by the volume number lookup equation
error_flag <- as.integer(0)  # for storing error messages
vol_equation <- as.character("")

vol_equation_function <- .Fortran("getvoleq_r",
                                  region,forest,district,species,
                                  vol_equation,error_flag)
vol_equation_function[[5]]


# 2 build on the above to obtain volume for tree

# define more inputs
vol_equation <- vol_equation_function[[5]]
PROD="01"            # 1-sawimber tree; 2-pulp; 3-roundwood
DBHOB=19.6
HTTOT=76
MTOPP=7             # Minimum top diameter inside bark for primary product
MTOPS=4             # Minimum top diameter inside bark for secondary product
STUMP=1             # stump ht
LIVE="L"            # live or dead
CTYPE="C"           # special cruise flag; C returns 0 if missing inputs
HT1PRD=0            # ht to min diameter for primary product
HT2PRD=0            # ht to min diameter for secondar product
UPSHT1=0            # Upper stem height in feet where upper stem diameter (0 if none)
UPSD1=0             # upper stem diameter at that ht (0 if none)
FCLASS=0            # girard form class
DBTBH=0             # double bark thickness at bh
BTR=0               # dib/dob *100
BRKHT=0             # ht to broken top
BRKHTD=0            # broken top diam
CR=0                # crown ratio
CULL=0              # pct cubic foot vol rotten or missing
DECAYCD=0           # decay code for standing dead
CULLMSTOP=0

# define the outputs
VOL=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
LOGVOL=matrix(0,7,20)
LOGDIA=matrix(0,21,3)
LOGLEN=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
BOLHT=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DRYBIO=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
GRNBIO=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
TLOGS=0
NOLOGP=0
NOLOGS=0
ERRFLG=0

# run the volume function
vol_func <- .Fortran("vollibnvb_r",
                     vol_equation,region,forest,district,species,
                     as.double(DBHOB),as.double(HTTOT),as.double(MTOPP),
                     as.double(MTOPS),as.double(HT1PRD),as.double(HT2PRD),as.double(UPSHT1),
                     as.double(UPSD1),as.double(STUMP),as.integer(FCLASS),as.double(DBTBH),as.double(BTR),
                     as.double(VOL),as.double(LOGVOL),as.double(LOGDIA),as.double(LOGLEN),as.double(BOLHT),
                     as.integer(TLOGS),as.double(NOLOGP),as.double(NOLOGS),error_flag,
                     as.double(BRKHT),as.double(BRKHTD),as.double(DRYBIO),as.double(GRNBIO),as.double(CR),
                     as.double(CULL),as.integer(DECAYCD), as.double(CULLMSTOP),as.character(CTYPE),
                     as.character(LIVE), as.character(PROD))
vol_func[[18]]

v.names <- c("total_cuft","gross_scrib","net_scrib",
             "gross_merch_cuft","net_merch_cuft","merch_cord",
             "gross_secnd_cuft","net_secnd_cuft","secnd_cord",
             "gross_intl","net_intl","gross_secnd_scrib","net_secnd_scrib",
             "stump_vol","tip_vol")
vols <- vol_func[[18]]; names(vols) <- v.names

vols[1]
vols[4]+vols[14]+vols[15] # not same as above because of merchandizing


# try a differnt 3-point volume equation
vol_equation
vol_equation2 <- "I00FW3W202"
UPSHT1=30
UPSD1=15

vol_func <- .Fortran("vollibnvb_r",
                     vol_equation2,region,forest,district,species,
                     as.double(DBHOB),as.double(HTTOT),as.double(MTOPP),
                     as.double(MTOPS),as.double(HT1PRD),as.double(HT2PRD),as.double(UPSHT1),
                     as.double(UPSD1),as.double(STUMP),as.integer(FCLASS),as.double(DBTBH),as.double(BTR),
                     as.double(VOL),as.double(LOGVOL),as.double(LOGDIA),as.double(LOGLEN),as.double(BOLHT),
                     as.integer(TLOGS),as.double(NOLOGP),as.double(NOLOGS),error_flag,
                     as.double(BRKHT),as.double(BRKHTD),as.double(DRYBIO),as.double(GRNBIO),as.double(CR),
                     as.double(CULL),as.integer(DECAYCD), as.double(CULLMSTOP),as.character(CTYPE),
                     as.character(LIVE), as.character(PROD))
vol_func[[18]][1] # new total cubic
vols[1]           # old total cubic


# 4 get taper profile

# dib at a given height 
taper <- data.frame(ht_above_ground = seq(1,floor(HTTOT),by=1),
                    dib_2pt = 0,
                    dib_3pt = 0)
for (i in 1:nrow(taper)){
  stemdib <- 0
  stemht <- taper$ht_above_ground[i]
  # 2 pt equation
  dibs <- .Fortran("calcdib_r",
                   vol_equation,region,forest,
                   as.double(DBHOB),as.double(HTTOT),
                   as.double(stemht),as.double(stemdib),
                   error_flag,as.double(UPSHT1),as.double(UPSD1))
  taper$dib_2pt[i] <- dibs[[7]]
  # 3 pt equation
  dibs <- .Fortran("calcdib_r",
                   vol_equation2,region,forest,
                   as.double(DBHOB),as.double(HTTOT),
                   as.double(stemht),as.double(stemdib),
                   error_flag,as.double(UPSHT1),as.double(UPSD1))
  taper$dib_3pt[i] <- dibs[[7]]
}

with(taper,plot(ht_above_ground,dib_2pt,type="l"))
with(taper,points(ht_above_ground,dib_3pt,type="l",col="red"))

