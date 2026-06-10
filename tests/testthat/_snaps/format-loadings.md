# format.LOADINGS aligns decimals and renders a plain table

    Code
      print(make_loadings())
    Output
                    F1     F2     F3 
      fun          .820  -.110   .050
      friends_lo   .450   .600  -.020
      enjoy       -.300   .710   .123
      hurt         .040  -.050  1.080
      

# format.LOADINGS prints communalities and a legend

    Code
      print(make_loadings(), h2 = h2, legend = TRUE)
    Output
                    F1     F2     F3     h2     u2 
      fun          .820  -.110   .050   .700   .300
      friends_lo   .450   .600  -.020   .580   .420
      enjoy       -.300   .710   .123   .630   .370
      hurt         .040  -.050  1.080  1.180  -.180
      
      Legend:
        bold = |loading| >= .300
        grey = below cutoff
        red h2/u2 = Heywood-relevant value
      

# format.LOADINGS sorts rows when requested

    Code
      print(make_loadings(), sort_loadings = "clustered")
    Output
                    F1     F2     F3 
      fun          .820  -.110   .050
      enjoy       -.300   .710   .123
      friends_lo   .450   .600  -.020
      hurt         .040  -.050  1.080
      

# wide matrices wrap into stacked blocks via max_factors_per_block

    Code
      print(make_wide_loadings(), max_factors_per_block = 3)
    Output
      Factors F1-F3 (block 1/3)
            F1     F2     F3 
      v1  -.900  -.670  -.440
      v2  -.840  -.610  -.380
      v3  -.780  -.550  -.320
      v4  -.730  -.490  -.260
      Factors F4-F6 (block 2/3)
            F4    F5    F6 
      v1  -.200  .030  .260
      v2  -.150  .090  .320
      v3  -.090  .150  .380
      v4  -.030  .200  .440
      Factors F7-F8 (block 3/3)
           F7    F8 
      v1  .490  .730
      v2  .550  .780
      v3  .610  .840
      v4  .670  .900
      

# wide matrices wrap into stacked blocks on a narrow console

    Code
      print(make_wide_loadings())
    Output
      Factors F1-F5 (block 1/2)
            F1     F2     F3     F4    F5 
      v1  -.900  -.670  -.440  -.200  .030
      v2  -.840  -.610  -.380  -.150  .090
      v3  -.780  -.550  -.320  -.090  .150
      v4  -.730  -.490  -.260  -.030  .200
      Factors F6-F8 (block 2/2)
           F6    F7    F8 
      v1  .260  .490  .730
      v2  .320  .550  .780
      v3  .380  .610  .840
      v4  .440  .670  .900
      

# format.SLLOADINGS flags a Heywood case

    Code
      print(make_sl(heywood = TRUE))
    Output
            g     F1     F2     h2     u2 
      i1  .700   .200   .100   .530   .470
      i2  .650  1.050  -.050  1.100  -.100
      i3  .400   .100   .550   .460   .540
      
      ! Results contain a Heywood case!
      

# format.SLLOADINGS prints cleanly without Heywood cases

    Code
      print(make_sl(heywood = FALSE))
    Output
            g    F1     F2    h2    u2 
      i1  .700  .200   .100  .530  .470
      i2  .650  .620  -.050  .700  .300
      i3  .400  .100   .550  .460  .540
      

