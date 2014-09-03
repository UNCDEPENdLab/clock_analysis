#!/usr/bin/env bash

#hard code header
hdr="lunaid AgeYearsDecimal sex scandate" # wverb_iq wperf_iq wfull4 wfull2 bdiTotal sssTAS sssES sssDIS sssB sssTOT asrTInternal asrTExternal asrTTotalProb"

echo $hdr | tr " " "\t" > subinfo_db

while read lunaid scandate; do

  #echo -ne "$lunaid	$age	$scandate	$rest	"
  mysql -h lncddb -u lncd lunadb_nightly -BNe "
select
  info.lunaid,
  TIMESTAMPDIFF(HOUR,info.DateOfBirth,lg.VisitDate)/8766.0 AS AgeYearsDecimal,
  case sexid when 1 then 'male' when 2 then 'female' else '?' end as sex,
  date_format(lg.VisitDate,'%Y-%m-%d') as scandate,
  dupps.*
from tSubjectInfo as info
left join tVisitLog as lg
   on lg.lunaid=info.lunaid
left join dupps
   on lg.lunaid=dupps.lunaid
where info.lunaid = '$lunaid'
  and date_format(lg.VisitDate,'%Y%m%d') = '$scandate'
order by  abs( datediff(dupps.visitdate,lg.visitdate) )
limit 1;
" 
echo
done < sublist_id_date | sed 's///g;/^$/d' >> subinfo_db

#  | sed  's///g;/^$/d' >> icaSubjInfo_DB.txt


#   dwasi.wverb_iq, dwasi.wperf_iq, dwasi.wfull4, dwasi.wfull2,
#   dBDI.bdiTotal,
#   dss.sssTAS, dss.sssES, dss.sssDIS, dss.sssBS, dss.sssTOT,
#   dasr.asrTInternal, dasr.asrTExternal, dasr.asrTTotalProb

# left join dSensationSeeking as dss
#    on dss.lunaID = lg.lunaID
# left join dwasi
#    on dwasi.VisitID = dss.VisitID
# left join dbdi
#    on dbdi.VisitID = dss.VisitID
# left join dasr
#    on dasr.VisitID = dss.VisitID
#order by abs( datediff(dss.visitdate,lg.visitdate) )
