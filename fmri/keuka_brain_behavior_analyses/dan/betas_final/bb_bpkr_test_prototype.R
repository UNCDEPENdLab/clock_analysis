library(parallel)
library(doParallel)
library(pbkrtest)

m_rewom_etheta <- lmer(value ~ rewFunc *  Stream * omission_early_theta  +  (1|id), df)
summary(m_rewom_etheta)
car::Anova(m_rewom_etheta, '3')
get_Lb_ddf(m_rewom_etheta, )
m_rewom_etheta_rs <- lmer(value ~ rewFunc *  Stream * omission_early_theta  +  (rewFunc|id), df)
summary(m_rewom_etheta_rs)
car::Anova(m_rewom_etheta_rs, '3')

m_rewom_etheta_rs0 <- lmer(value ~ rewFunc *  Stream +  (rewFunc|id), df)
summary(m_rewom_etheta_rs0)
car::Anova(m_rewom_etheta_rs0, '3')




KRmodcomp(m_rewom_etheta_rs, m_rewom_etheta_rs0)
getKR(
  object,
  name = c("ndf", "ddf", "Fstat", "p.value", "F.scaling", "FstatU", "p.valueU", "aux")
)



f <- Sys.getenv('PBS_NODEFILE')

ncores <- detectCores()
nodelist <- if (nzchar(f)) readLines(f) else rep('localhost', ncores)

cat("Node list allocated to this job\n")
print(nodelist)

cl <- makePSOCKcluster(nodelist, outfile='')
print(cl) ##; print(unclass(cl))
registerDoParallel(cl)

pb <- pbkrtest::PBmodcomp(m_rewom_etheta_rs, m_rewom_etheta_rs0, cl = cl)

getKR(
  pb,
  name = c("ndf", "ddf", "Fstat", "p.value", "F.scaling", "FstatU", "p.valueU", "aux")
)
