function h = eurCalendar
% Default target non-trading day data from 2000 to 2065
% 
% RondPoint 2012
% Last Modified: 21.03.2012 R. Baviera, A. Cassaro 

h = [...
    730486 ;... % 1-Jan-2000
    730597 ;... % 21-Apr-2000
    730600 ;... % 24-Apr-2000
    730607 ;... % 1-May-2000
    730845 ;... % 25-Dec-2000
    730846 ;... % 26-Dec-2000
    730852 ;... % 1-Jan-2001
    730954 ;... % 13-Apr-2001
    730957 ;... % 16-Apr-2001
    730972 ;... % 1-May-2001
    731210 ;... % 25-Dec-2001
    731211 ;... % 26-Dec-2001
    731217 ;... % 1-Jan-2002
    731304 ;... % 29-Mar-2002
    731307 ;... % 1-Apr-2002
    731337 ;... % 1-May-2002
    731575 ;... % 25-Dec-2002
    731576 ;... % 26-Dec-2002
    731582 ;... % 1-Jan-2003
    731689 ;... % 18-Apr-2003
    731692 ;... % 21-Apr-2003
    731702 ;... % 1-May-2003
    731940 ;... % 25-Dec-2003
    731941 ;... % 26-Dec-2003
    731947 ;... % 1-Jan-2004
    732046 ;... % 9-Apr-2004
    732049 ;... % 12-Apr-2004
    732068 ;... % 1-May-2004
    732306 ;... % 25-Dec-2004
    732307 ;... % 26-Dec-2004
    732313 ;... % 1-Jan-2005
    732396 ;... % 25-Mar-2005
    732399 ;... % 28-Mar-2005
    732433 ;... % 1-May-2005
    732671 ;... % 25-Dec-2005
    732672 ;... % 26-Dec-2005
    732678 ;... % 1-Jan-2006
    732781 ;... % 14-Apr-2006
    732784 ;... % 17-Apr-2006
    732798 ;... % 1-May-2006
    733036 ;... % 25-Dec-2006
    733037 ;... % 26-Dec-2006
    733043 ;... % 1-Jan-2007
    733138 ;... % 6-Apr-2007
    733141 ;... % 9-Apr-2007
    733163 ;... % 1-May-2007
    733401 ;... % 25-Dec-2007
    733402 ;... % 26-Dec-2007
    733408 ;... % 1-Jan-2008
    733488 ;... % 21-Mar-2008
    733491 ;... % 24-Mar-2008
    733529 ;... % 1-May-2008
    733767 ;... % 25-Dec-2008
    733768 ;... % 26-Dec-2008
    733774 ;... % 1-Jan-2009
    733873 ;... % 10-Apr-2009
    733876 ;... % 13-Apr-2009
    733894 ;... % 1-May-2009
    734132 ;... % 25-Dec-2009
    734133 ;... % 26-Dec-2009
    734139 ;... % 1-Jan-2010
    734230 ;... % 2-Apr-2010
    734233 ;... % 5-Apr-2010
    734259 ;... % 1-May-2010
    734497 ;... % 25-Dec-2010
    734498 ;... % 26-Dec-2010
    734504 ;... % 1-Jan-2011
    734615 ;... % 22-Apr-2011
    734618 ;... % 25-Apr-2011
    734624 ;... % 1-May-2011
    734862 ;... % 25-Dec-2011
    734863 ;... % 26-Dec-2011
    734869 ;... % 1-Jan-2012
    734965 ;... % 6-Apr-2012
    734968 ;... % 9-Apr-2012
    734990 ;... % 1-May-2012
    735228 ;... % 25-Dec-2012
    735229 ;... % 26-Dec-2012
    735235 ;... % 1-Jan-2013
    735322 ;... % 29-Mar-2013
    735325 ;... % 1-Apr-2013
    735355 ;... % 1-May-2013
    735593 ;... % 25-Dec-2013
    735594 ;... % 26-Dec-2013
    735600 ;... % 1-Jan-2014
    735707 ;... % 18-Apr-2014
    735710 ;... % 21-Apr-2014
    735720 ;... % 1-May-2014
    735958 ;... % 25-Dec-2014
    735959 ;... % 26-Dec-2014
    735965 ;... % 1-Jan-2015
    736057 ;... % 3-Apr-2015
    736060 ;... % 6-Apr-2015
    736085 ;... % 1-May-2015
    736323 ;... % 25-Dec-2015
    736324 ;... % 26-Dec-2015
    736330 ;... % 1-Jan-2016
    736414 ;... % 25-Mar-2016
    736417 ;... % 28-Mar-2016
    736451 ;... % 1-May-2016
    736689 ;... % 25-Dec-2016
    736690 ;... % 26-Dec-2016
    736696 ;... % 1-Jan-2017
    736799 ;... % 14-Apr-2017
    736802 ;... % 17-Apr-2017
    736816 ;... % 1-May-2017
    737054 ;... % 25-Dec-2017
    737055 ;... % 26-Dec-2017
    737061 ;... % 1-Jan-2018
    737149 ;... % 30-Mar-2018
    737152 ;... % 2-Apr-2018
    737181 ;... % 1-May-2018
    737419 ;... % 25-Dec-2018
    737420 ;... % 26-Dec-2018
    737426 ;... % 1-Jan-2019
    737534 ;... % 19-Apr-2019
    737537 ;... % 22-Apr-2019
    737546 ;... % 1-May-2019
    737784 ;... % 25-Dec-2019
    737785 ;... % 26-Dec-2019
    737791 ;... % 1-Jan-2020
    737891 ;... % 10-Apr-2020
    737894 ;... % 13-Apr-2020
    737912 ;... % 1-May-2020
    738150 ;... % 25-Dec-2020
    738151 ;... % 26-Dec-2020
    738157 ;... % 1-Jan-2021
    738248 ;... % 2-Apr-2021
    738251 ;... % 5-Apr-2021
    738277 ;... % 1-May-2021
    738515 ;... % 25-Dec-2021
    738516 ;... % 26-Dec-2021
    738522 ;... % 1-Jan-2022
    738626 ;... % 15-Apr-2022
    738629 ;... % 18-Apr-2022
    738642 ;... % 1-May-2022
    738880 ;... % 25-Dec-2022
    738881 ;... % 26-Dec-2022
    738887 ;... % 1-Jan-2023
    738983 ;... % 7-Apr-2023
    738986 ;... % 10-Apr-2023
    739007 ;... % 1-May-2023
    739245 ;... % 25-Dec-2023
    739246 ;... % 26-Dec-2023
    739252 ;... % 1-Jan-2024
    739340 ;... % 29-Mar-2024
    739343 ;... % 1-Apr-2024
    739373 ;... % 1-May-2024
    739611 ;... % 25-Dec-2024
    739612 ;... % 26-Dec-2024
    739618 ;... % 1-Jan-2025
    739725 ;... % 18-Apr-2025
    739728 ;... % 21-Apr-2025
    739738 ;... % 1-May-2025
    739976 ;... % 25-Dec-2025
    739977 ;... % 26-Dec-2025
    739983 ;... % 1-Jan-2026
    740075 ;... % 3-Apr-2026
    740078 ;... % 6-Apr-2026
    740103 ;... % 1-May-2026
    740341 ;... % 25-Dec-2026
    740342 ;... % 26-Dec-2026
    740348 ;... % 1-Jan-2027
    740432 ;... % 26-Mar-2027
    740435 ;... % 29-Mar-2027
    740468 ;... % 1-May-2027
    740706 ;... % 25-Dec-2027
    740707 ;... % 26-Dec-2027
    740713 ;... % 1-Jan-2028
    740817 ;... % 14-Apr-2028
    740820 ;... % 17-Apr-2028
    740834 ;... % 1-May-2028
    741072 ;... % 25-Dec-2028
    741073 ;... % 26-Dec-2028
    741079 ;... % 1-Jan-2029
    741167 ;... % 30-Mar-2029
    741170 ;... % 2-Apr-2029
    741199 ;... % 1-May-2029
    741437 ;... % 25-Dec-2029
    741438 ;... % 26-Dec-2029
    741444 ;... % 1-Jan-2030
    741552 ;... % 19-Apr-2030
    741555 ;... % 22-Apr-2030
    741564 ;... % 1-May-2030
    741802 ;... % 25-Dec-2030
    741803 ;... % 26-Dec-2030
    741809 ;... % 1-Jan-2031
    741909 ;... % 11-Apr-2031
    741912 ;... % 14-Apr-2031
    741929 ;... % 1-May-2031
    742167 ;... % 25-Dec-2031
    742168 ;... % 26-Dec-2031
    742174 ;... % 1-Jan-2032
    742259 ;... % 26-Mar-2032
    742262 ;... % 29-Mar-2032
    742295 ;... % 1-May-2032
    742533 ;... % 25-Dec-2032
    742534 ;... % 26-Dec-2032
    742540 ;... % 1-Jan-2033
    742644 ;... % 15-Apr-2033
    742647 ;... % 18-Apr-2033
    742660 ;... % 1-May-2033
    742898 ;... % 25-Dec-2033
    742899 ;... % 26-Dec-2033
    742905 ;... % 1-Jan-2034
    743001 ;... % 7-Apr-2034
    743004 ;... % 10-Apr-2034
    743025 ;... % 1-May-2034
    743263 ;... % 25-Dec-2034
    743264 ;... % 26-Dec-2034
    743270 ;... % 1-Jan-2035
    743351 ;... % 23-Mar-2035
    743354 ;... % 26-Mar-2035
    743390 ;... % 1-May-2035
    743628 ;... % 25-Dec-2035
    743629 ;... % 26-Dec-2035
    743635 ;... % 1-Jan-2036
    743736 ;... % 11-Apr-2036
    743739 ;... % 14-Apr-2036
    743756 ;... % 1-May-2036
    743994 ;... % 25-Dec-2036
    743995 ;... % 26-Dec-2036
    744001 ;... % 1-Jan-2037
    744093 ;... % 3-Apr-2037
    744096 ;... % 6-Apr-2037
    744121 ;... % 1-May-2037
    744359 ;... % 25-Dec-2037
    744360 ;... % 26-Dec-2037
    744366 ;... % 1-Jan-2038
    744478 ;... % 23-Apr-2038
    744481 ;... % 26-Apr-2038
    744486 ;... % 1-May-2038
    744724 ;... % 25-Dec-2038
    744725 ;... % 26-Dec-2038
    744731 ;... % 1-Jan-2039
    744828 ;... % 8-Apr-2039
    744831 ;... % 11-Apr-2039
    744851 ;... % 1-May-2039
    745089 ;... % 25-Dec-2039
    745090 ;... % 26-Dec-2039
    745096 ;... % 1-Jan-2040
    745185 ;... % 30-Mar-2040
    745188 ;... % 2-Apr-2040
    745217 ;... % 1-May-2040
    745455 ;... % 25-Dec-2040
    745456 ;... % 26-Dec-2040
    745462 ;... % 1-Jan-2041
    745570 ;... % 19-Apr-2041
    745573 ;... % 22-Apr-2041
    745582 ;... % 1-May-2041
    745820 ;... % 25-Dec-2041
    745821 ;... % 26-Dec-2041
    745827 ;... % 1-Jan-2042
    745920 ;... % 4-Apr-2042
    745923 ;... % 7-Apr-2042
    745947 ;... % 1-May-2042
    746185 ;... % 25-Dec-2042
    746186 ;... % 26-Dec-2042
    746192 ;... % 1-Jan-2043
    746277 ;... % 27-Mar-2043
    746280 ;... % 30-Mar-2043
    746312 ;... % 1-May-2043
    746550 ;... % 25-Dec-2043
    746551 ;... % 26-Dec-2043
    746557 ;... % 1-Jan-2044
    746662 ;... % 15-Apr-2044
    746665 ;... % 18-Apr-2044
    746678 ;... % 1-May-2044
    746916 ;... % 25-Dec-2044
    746917 ;... % 26-Dec-2044
    746923 ;... % 1-Jan-2045
    747019 ;... % 7-Apr-2045
    747022 ;... % 10-Apr-2045
    747043 ;... % 1-May-2045
    747281 ;... % 25-Dec-2045
    747282 ;... % 26-Dec-2045
    747288 ;... % 1-Jan-2046
    747369 ;... % 23-Mar-2046
    747372 ;... % 26-Mar-2046
    747408 ;... % 1-May-2046
    747646 ;... % 25-Dec-2046
    747647 ;... % 26-Dec-2046
    747653 ;... % 1-Jan-2047
    747754 ;... % 12-Apr-2047
    747757 ;... % 15-Apr-2047
    747773 ;... % 1-May-2047
    748011 ;... % 25-Dec-2047
    748012 ;... % 26-Dec-2047
    748018 ;... % 1-Jan-2048
    748111 ;... % 3-Apr-2048
    748114 ;... % 6-Apr-2048
    748139 ;... % 1-May-2048
    748377 ;... % 25-Dec-2048
    748378 ;... % 26-Dec-2048
    748384 ;... % 1-Jan-2049
    748489 ;... % 16-Apr-2049
    748492 ;... % 19-Apr-2049
    748504 ;... % 1-May-2049
    748742 ;... % 25-Dec-2049
    748743 ;... % 26-Dec-2049
    748749 ;... % 1-Jan-2050
    748846 ;... % 8-Apr-2050
    748849 ;... % 11-Apr-2050
    748869 ;... % 1-May-2050
    749107 ;... % 25-Dec-2050
    749108 ;... % 26-Dec-2050
    749114 ;... % 1-Jan-2051
    749203 ;... % 31-Mar-2051
    749206 ;... % 3-Apr-2051
    749234 ;... % 1-May-2051
    749472 ;... % 25-Dec-2051
    749473 ;... % 26-Dec-2051
    749479 ;... % 1-Jan-2052
    749588 ;... % 19-Apr-2052
    749591 ;... % 22-Apr-2052
    749600 ;... % 1-May-2052
    749838 ;... % 25-Dec-2052
    749839 ;... % 26-Dec-2052
    749845 ;... % 1-Jan-2053
    749938 ;... % 4-Apr-2053
    749941 ;... % 7-Apr-2053
    749965 ;... % 1-May-2053
    750203 ;... % 25-Dec-2053
    750204 ;... % 26-Dec-2053
    750210 ;... % 1-Jan-2054
    750295 ;... % 27-Mar-2054
    750298 ;... % 30-Mar-2054
    750330 ;... % 1-May-2054
    750568 ;... % 25-Dec-2054
    750569 ;... % 26-Dec-2054
    750575 ;... % 1-Jan-2055
    750680 ;... % 16-Apr-2055
    750683 ;... % 19-Apr-2055
    750695 ;... % 1-May-2055
    750933 ;... % 25-Dec-2055
    750934 ;... % 26-Dec-2055
    750940 ;... % 1-Jan-2056
    751030 ;... % 31-Mar-2056
    751033 ;... % 3-Apr-2056
    751061 ;... % 1-May-2056
    751299 ;... % 25-Dec-2056
    751300 ;... % 26-Dec-2056
    751306 ;... % 1-Jan-2057
    751415 ;... % 20-Apr-2057
    751418 ;... % 23-Apr-2057
    751426 ;... % 1-May-2057
    751664 ;... % 25-Dec-2057
    751665 ;... % 26-Dec-2057
    751671 ;... % 1-Jan-2058
    751772 ;... % 12-Apr-2058
    751775 ;... % 15-Apr-2058
    751791 ;... % 1-May-2058
    752029 ;... % 25-Dec-2058
    752030 ;... % 26-Dec-2058
    752036 ;... % 1-Jan-2059
    752122 ;... % 28-Mar-2059
    752125 ;... % 31-Mar-2059
    752156 ;... % 1-May-2059
    752394 ;... % 25-Dec-2059
    752395 ;... % 26-Dec-2059
    752401 ;... % 1-Jan-2060
    752507 ;... % 16-Apr-2060
    752510 ;... % 19-Apr-2060
    752522 ;... % 1-May-2060
    752760 ;... % 25-Dec-2060
    752761 ;... % 26-Dec-2060
    752767 ;... % 1-Jan-2061
    752864 ;... % 8-Apr-2061
    752867 ;... % 11-Apr-2061
    752887 ;... % 1-May-2061
    753125 ;... % 25-Dec-2061
    753126 ;... % 26-Dec-2061
    753132 ;... % 1-Jan-2062
    753214 ;... % 24-Mar-2062
    753217 ;... % 27-Mar-2062
    753252 ;... % 1-May-2062
    753490 ;... % 25-Dec-2062
    753491 ;... % 26-Dec-2062
    753497 ;... % 1-Jan-2063
    753599 ;... % 13-Apr-2063
    753602 ;... % 16-Apr-2063
    753617 ;... % 1-May-2063
    753855 ;... % 25-Dec-2063
    753856 ;... % 26-Dec-2063
    753862 ;... % 1-Jan-2064
    753956 ;... % 4-Apr-2064
    753959 ;... % 7-Apr-2064
    753983 ;... % 1-May-2064
    754221 ;... % 25-Dec-2064
    754222 ;... % 26-Dec-2064
    754228 ;... % 1-Jan-2065
    754313 ;... % 27-Mar-2065
    754316 ;... % 30-Mar-2065
    754348 ;... % 1-May-2065
    754586 ;... % 25-Dec-2065
    754587 ;... % 26-Dec-2065
    ];

end % function