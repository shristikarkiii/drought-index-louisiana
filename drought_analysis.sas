/*---------------------------------------------------------------------------
  Louisiana Drought Index Analysis -- SAS
  Shristi Karki

  Reads the gap-free monthly precipitation grid produced by compute_spi.py,
  recomputes the 12-month SPI independently in SAS, cross-checks it against the
  Python result, and produces the station and seasonal summaries and plots.

  Recomputing rather than re-reading is deliberate: it makes this an
  independent check of the Python pipeline, not a restatement of it.

  Requires Base SAS only. Set PROJ below to the project folder.
---------------------------------------------------------------------------*/

%let proj = C:\Users\shris_gq83\Desktop\drought-louisiana;
libname dr "&proj.\data";

/*-- 1. import the monthly grid (station, year, month, prcp_mm, spi3/6/12) --*/
proc import datafile="&proj.\data\spi.csv" out=work.spi_in dbms=csv replace;
    getnames=yes; guessingrows=max;
run;

proc sort data=work.spi_in out=work.grid;
    by station_id year month;
run;

/*-- 2. 12-month accumulation. A window containing any missing month stays
       missing, so a gap is never silently treated as a dry month.          --*/
data work.acc;
    set work.grid;
    by station_id;
    array hist{12} _temporary_;
    retain nobs 0;
    if first.station_id then do;
        nobs = 0;
        do i = 1 to 12; hist{i} = .; end;
    end;
    nobs + 1;
    hist{ mod(nobs-1, 12) + 1 } = prcp_mm;
    if nobs >= 12 and cmiss(of hist{*}) = 0 then acc12 = sum(of hist{*});
    else acc12 = .;
    drop i;
run;

/*-- 3. gamma parameters per station and calendar month, fitted on the
       non-zero accumulations only (zeros handled separately below)        --*/
data work.pos work.zeros;
    set work.acc;
    if acc12 > . then do;
        if acc12 > 0 then output work.pos;
        else output work.zeros;
    end;
run;

proc sort data=work.pos; by station_id month; run;

proc univariate data=work.pos noprint;
    by station_id month;
    var acc12;
    histogram acc12 / gamma(theta=0 fitmethod=mle) outfit=work.gfit;
run;

/* probability of a zero accumulation, needed for the mixed distribution */
proc sql;
    create table work.qzero as
    select station_id, month,
           sum(case when acc12 = 0 then 1 else 0 end) / count(*) as q_zero
    from work.acc
    where acc12 > .
    group by station_id, month;
quit;

data work.params;
    merge work.gfit(keep=station_id month _SHAPE1_ _SCALE_
                    rename=(_SHAPE1_=shape _SCALE_=scale))
          work.qzero;
    by station_id month;
run;

/*-- 4. SPI: map the mixed CDF onto the standard normal --------------------*/
proc sort data=work.acc; by station_id month; run;

data work.spi_sas;
    merge work.acc work.params;
    by station_id month;
    if acc12 > . and shape > . and scale > . then do;
        cdf_raw = q_zero + (1 - q_zero) * cdf('GAMMA', acc12, shape, scale);
        cdf_raw = min(max(cdf_raw, 1e-6), 1 - 1e-6);   /* keep probit finite */
        spi12_sas = probit(cdf_raw);
    end;

    length drought_class $ 8;
    if      spi12_sas = .        then drought_class = '';
    else if spi12_sas <= -2.0    then drought_class = 'Extreme';
    else if spi12_sas <= -1.5    then drought_class = 'Severe';
    else if spi12_sas <= -1.0    then drought_class = 'Moderate';
    else                              drought_class = 'None';

    length season $ 12;
    select (month);
        when (12, 1, 2) season = 'Winter';
        when (3, 4, 5)  season = 'Spring';
        when (6, 7, 8)  season = 'Summer';
        otherwise       season = 'Autumn';
    end;

    diff = spi12_sas - spi12;    /* cross-check against the Python result */
run;

/*-- 5. cross-check ---------------------------------------------------------*/
title "Agreement between the SAS and Python SPI-12";
proc means data=work.spi_sas n mean std min max maxdec=4;
    var diff;
run;

proc corr data=work.spi_sas pearson;
    var spi12 spi12_sas;
run;
title;

/*-- 6. summaries -----------------------------------------------------------*/
title "Drought frequency, duration and severity by station (SPI-12)";
proc means data=work.spi_sas n mean std min p10 median maxdec=2;
    class station;
    var spi12_sas;
run;

proc freq data=work.spi_sas;
    tables station * drought_class / nocum nopercent norow;
    where drought_class ne '';
run;

title "Drought incidence by season";
proc freq data=work.spi_sas;
    tables season * drought_class / nocum nopercent norow;
    where drought_class ne '';
run;
title;

/*-- 7. drought events: runs of consecutive months at or below -1 -----------*/
proc sort data=work.spi_sas out=work.seq; by station_id year month; run;

data work.events;
    set work.seq;
    by station_id;
    retain run_len 0 sev 0 peak . start_y . start_m .;
    in_drought = (spi12_sas ne . and spi12_sas <= -1.0);

    if first.station_id then do; run_len = 0; sev = 0; peak = .; end;

    if in_drought then do;
        if run_len = 0 then do; start_y = year; start_m = month; peak = spi12_sas; end;
        run_len + 1;
        sev = sev - spi12_sas;                 /* accumulated deficit */
        peak = min(peak, spi12_sas);
    end;
    else if run_len > 0 then do;
        duration_months = run_len;
        severity = sev;
        peak_spi = peak;
        output;
        run_len = 0; sev = 0; peak = .;
    end;

    if last.station_id and run_len > 0 then do;
        duration_months = run_len; severity = sev; peak_spi = peak; output;
        run_len = 0; sev = 0; peak = .;
    end;

    keep station station_id start_y start_m duration_months severity peak_spi;
run;

title "Most severe drought episodes on record";
proc sort data=work.events; by descending severity; run;
proc print data=work.events(obs=10) noobs;
    var station start_y start_m duration_months severity peak_spi;
    format severity peak_spi 8.2;
run;
title;

/*-- 8. plots ---------------------------------------------------------------*/
ods graphics on / width=10in height=4in;

title "Statewide mean 12-month SPI";
proc sql;
    create table work.statewide as
    select year + (month - 1) / 12 as t, mean(spi12_sas) as mean_spi
    from work.spi_sas where spi12_sas > . group by t;
quit;

proc sgplot data=work.statewide noautolegend;
    band x=t lower=0 upper=mean_spi / fillattrs=(color=cx8c3b1e) transparency=.15;
    series x=t y=mean_spi / lineattrs=(color=cx3f3f3f thickness=1);
    refline -1 / axis=y lineattrs=(pattern=shortdash color=cx333333)
                 label="drought threshold";
    xaxis label="Year"; yaxis label="Mean SPI-12";
run;

title "Percent of months in drought by station";
proc sgplot data=work.spi_sas;
    hbar station / response=in_drought stat=mean
          fillattrs=(color=cx8c3b1e);
    refline 0.1587 / axis=x lineattrs=(pattern=shortdash) label="expected 15.9%";
    xaxis label="Share of months with SPI-12 <= -1";
    yaxis label=" " discreteorder=data;
run;

title "SPI-12 by season";
proc sgplot data=work.spi_sas;
    vbox spi12_sas / category=season fillattrs=(color=cxe2a76b);
    refline -1 / axis=y lineattrs=(pattern=shortdash) label="drought threshold";
    yaxis label="SPI-12"; xaxis label="Season";
run;
title;

ods graphics off;

/*-- 9. save the SAS result for comparison ---------------------------------*/
proc export data=work.spi_sas(keep=station_id station year month prcp_mm
                                   acc12 spi12 spi12_sas diff drought_class season)
    outfile="&proj.\data\spi_sas_check.csv" dbms=csv replace;
run;
