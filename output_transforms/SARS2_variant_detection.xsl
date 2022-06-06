<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xs="http://www.w3.org/2001/XMLSchema" xmlns:exsl="http://exslt.org/common" xmlns:str="http://exslt.org/strings" extension-element-prefixes="exsl str">
    <xsl:import href="http://exslt.org/str/functions/replace/str.replace.function.xsl"/>
    <xsl:output method="xhtml" doctype-public="-//W3C//DTD XHTML 1.0 Strict//EN" doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd" omit-xml-declaration="yes" encoding="UTF-8" indent="yes"/>
    <xsl:template match="/analysis">
        <html>
        <head>
            <title>Run Summary for: <xsl:value-of select="@run_name"/></title>
            <style>
                .freeze-table {
                    border-spacing: 0;
                    font-size: 14px;
                    padding: 0;
                    border: 1px solid #ccc;
                }
                th {
                    top: 0;
                    position: sticky;
                    background-color: #666;
                    color: #fff;
                    z-index: 20;
                    min-height: 30px;
                    height: 30px;
                    text-align: center;
                    border: 1px solid #fff;
                }
                .second-header th{
                    top: 30px;
                }
                td {
                    border: 1px solid #ccc;
                    background-color: #fff;
                }
                .fixed-header {
                    z-index: 50;
                }
                .fixed-col-sample {
                    left: 0;
                    position: sticky;
                    width: 240px;
                    border-right: 3px solid #ccc
                }
            </style>
        </head>
        <body>
            <center><h1>ASAP Run Summary for: <xsl:value-of select="@run_name"/></h1></center>
	    <br/>
            <em>Run stats</em>
            <table>
                <tr>
                <td>Total Samples</td>
                <td><xsl:value-of select="count(//sample[not(contains(@name, 'NTC') or contains(@name, 'NEC') or contains(@name, 'PTC'))])"/></td>
                </tr>
                <tr>
                <td>Total (&gt;90%)</td>
                <td><xsl:value-of select="count(//sample[not(contains(@name, 'NTC') or contains(@name, 'NEC') or contains(@name, 'PTC'))]//breadth[text() &gt;= 90])"/></td>
                </tr>
                <tr>
                <td>Total (80-90%)</td>
                <td><xsl:value-of select="count(//sample[not(contains(@name, 'NTC') or contains(@name, 'NEC') or contains(@name, 'PTC'))]//breadth[text() &gt;= 80 and text() &lt; 90])"/></td>
                </tr>
                <tr>
                <td>Total (&lt;80%)</td>
                <td><xsl:value-of select="count(//sample[not(contains(@name, 'NTC') or contains(@name, 'NEC') or contains(@name, 'PTC'))]//breadth[text() &lt; 80])"/></td>
                </tr>
            </table>
            <br/>
	    <em>Alignment stats for each sample</em>
            <table class="freeze-table" width="100%">
                <thead>
	        <tr>
	    	<th class="fixed-col-sample fixed-header">Sample</th>
                <th>Total Reads</th>
	    	<xsl:for-each select="sample[1]">
	    	    <xsl:for-each select="assay">
	    	        <th><xsl:value-of select='@name'/> Aligned</th>
                        <th>% Aligned</th>
                        <th>Coverage Breadth</th>
                        <th>Average Depth</th>
                        <th># SNPs + INDELs</th>
	    	    </xsl:for-each>
	    	</xsl:for-each>
	    	</tr>
                </thead>
                <xsl:for-each select="sample">
                    <tr>
                    <td nowrap="true" class="fixed-col-sample"><a href="{/analysis/@run_name}/{./@name}.html"><xsl:value-of select="@name"/></a></td>
                    <xsl:variable name="TOTAL_READS" select="@mapped_reads + @unmapped_reads"/>
                    <td align="right"><xsl:value-of select="$TOTAL_READS"/></td>
		    <xsl:for-each select="assay">
		        <td align="right"><xsl:value-of select="amplicon/@reads"/></td>
	    		<td align="right"><xsl:value-of select='format-number(amplicon/@reads div $TOTAL_READS * 100, "##.##")'/>%</td>
	    		<td align="right"><xsl:value-of select='format-number(amplicon/breadth, "##.##")'/>%</td>
	    		<td align="right"><xsl:value-of select='format-number(amplicon/average_depth, "##.##")'/></td>
                        <td align="right"><xsl:value-of select='count(amplicon/snp/snp_call[@percent &gt;= 10])'/></td>
		    </xsl:for-each>
                    </tr>
                    <xsl:apply-templates select="."/>
                </xsl:for-each>
            </table>
            <br/>
            <a href="{@run_name}_variants.html">Click here for the variant report</a>
        </body>
        </html>

<!-- Variant Report -->
    <exsl:document method="html" href="{@run_name}_variants.html">
        <html>
        <head>
            <title>Variant Report for: <xsl:value-of select="@run_name"/></title>
            <style>
                .freeze-table {
                    border-spacing: 0;
                    font-size: 14px;
                    padding: 0;
                    border: 1px solid #ccc;
                }
                th {
                    top: 0;
                    position: sticky;
                    background-color: #666;
                    color: #fff;
                    z-index: 20;
                    min-height: 30px;
                    height: 30px;
                    text-align: center;
                    border: 1px solid #fff;
                }
                .second-header th{
                    top: 30px;
                }
                td {
                    border: 1px solid #ccc;
                    background-color: #fff;
                }
                .fixed-header {
                    z-index: 50;
                }
                .fixed-col-sample {
                    left: 0;
                    position: sticky;
                    width: 240px;
                }
                .fixed-col-breadth {
                    left: 240px;
                    position: sticky;
                    border-right: 3px solid #ccc
                }
            </style>
        </head>
        <body>
            <center><h1>COVID Variant Report for: <xsl:value-of select="@run_name"/></h1></center>
	    <br/>
            <table class="freeze-table" width="100%">
                <thead>
	        <tr>
	    	<th rowspan="2" class="fixed-col-sample fixed-header">Sample (Total: <xsl:value-of select="count(//sample[not(contains(@name, 'NTC'))])"/>)</th>
	    	<th rowspan="2" class="fixed-col-breadth fixed-header">Coverage Breadth</th>
                <th colspan="13">Delta Variants</th>
                <th colspan="1">AY.3</th>
                <th colspan="2">AY.4</th>
                <th colspan="1">AY.20</th>
                <th colspan="1">AY.24</th>
                <th colspan="2">AY.25</th>
                <th colspan="3">AY.26</th>
                <th colspan="2">AY.39</th>
                <th colspan="2">AY.44</th>
                <th colspan="2">AY.103</th>
                <th colspan="5">Gamma Variant</th>
                <th colspan="4">Alpha Variant</th>
                <th colspan="4">Beta Variant</th>
                <th colspan="4">Other Spike Mutations</th>
                <th colspan="2">Not A.1</th>
                <th colspan="2">A.1</th>
	    	<th rowspan="2">All Deletions</th>
                </tr>
                <tr class="second-header">
                <th>S:T19R</th>
                <th>S:G142D</th>
                <th>S:L452R</th>
                <th>S:T478K</th>
                <th>S:P681R</th>
                <th>S:D950N</th>
                <th>ORF3a:S26L</th>
                <th>M:I82T</th>
                <th>ORF7a:V82A</th>
                <th>ORF7a:T120I</th>
                <th>N:D63G</th>
                <th>N:R203M</th>
                <th>N:D377Y</th>
                <th>ORF1a:I3731V</th>
                <th>ORF1a:A2529V</th>
                <th>S:T95I</th>
                <th>S:V1104L</th>
                <th>S:A222V</th>
                <th>ORF3a:E239Q</th>
                <th>ORF7a:G38syn</th>
                <th>S:A222V</th>
                <th>S:V1264L</th>
                <th>ORF6:K48N</th>
                <th>ORF1b:Q2635H</th>
                <th>ORF7a:V71I</th>
                <th>ORF1a:F1925syn</th>
                <th>ORF1a:A2554V</th>
                <th>ORF1b:I1257V</th>
                <th>S:I882syn</th>
                <th>S:L18F</th>
                <th>S:K417T</th>
                <th>S:E484K</th>
                <th>S:N501Y</th>
                <th>S:H655Y</th>
                <th>S:N501Y</th>
                <th>S:A570D</th>
                <th>S:S982A</th>
                <th>S:D1118H</th>
                <th>S:K417N</th>
                <th>S:E484K</th>
                <th>S:N501Y</th>
                <th>S:A701V</th>
                <th>S:K378E</th>
                <th>S:K378N</th>
                <th>S:G446D</th>
                <th>S:Q498R</th>
                <th>S:D614G</th>
                <th>ORF1b:P314L</th>
                <th>ORF1b:P1427L</th>
                <th>ORF1b:Y1464C</th>
                </tr>
                </thead>
                <xsl:for-each select="sample">
                    <xsl:variable name="TOTAL_READS" select="@mapped_reads + @unmapped_reads"/>
                    <xsl:variable name="BG_T19R"><xsl:choose><xsl:when test="descendant::snp[@name='S:T19R']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:T19R']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:T19R']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:T19R']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_G142D"><xsl:choose><xsl:when test="descendant::snp[@name='S:G142D']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:G142D']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:G142D']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:G142D']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_L452R"><xsl:choose><xsl:when test="descendant::snp[@name='S:L452R']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:L452R']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:L452R']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:L452R']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_T478K"><xsl:choose><xsl:when test="descendant::snp[@name='S:T478K']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:T478K']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:T478K']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:T478K']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_P681R"><xsl:choose><xsl:when test="descendant::snp[@name='S:P681R']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:P681R']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:P681R']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:P681R']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_D950N"><xsl:choose><xsl:when test="descendant::snp[@name='S:D950N']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:D950N']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:D950N']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:D950N']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_S26L"><xsl:choose><xsl:when test="descendant::snp[@name='ORF3a:S26L']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF3a:S26L']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='ORF3a:S26L']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF3a:S26L']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_I82T"><xsl:choose><xsl:when test="descendant::snp[@name='M:I82T']/snp_call/@percent &gt; 80 and descendant::snp[@name='M:I82T']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='M:I82T']/snp_call/@percent &gt; 80 and descendant::snp[@name='M:I82T']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_V82A"><xsl:choose><xsl:when test="descendant::snp[@name='ORF7a:V82A']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF7a:V82A']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='ORF7a:V82A']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF7a:V82A']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_T120I"><xsl:choose><xsl:when test="descendant::snp[@name='ORF7a:T120I']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF7a:T120I']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='ORF7a:T120I']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF7a:T120I']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_D63G"><xsl:choose><xsl:when test="descendant::snp[@name='N:D63G']/snp_call/@percent &gt; 80 and descendant::snp[@name='N:D63G']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='N:D63G']/snp_call/@percent &gt; 80 and descendant::snp[@name='N:D63G']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_R203M"><xsl:choose><xsl:when test="descendant::snp[@name='N:R203M']/snp_call/@percent &gt; 80 and descendant::snp[@name='N:R203M']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='N:R203M']/snp_call/@percent &gt; 80 and descendant::snp[@name='N:R203M']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_D377Y"><xsl:choose><xsl:when test="descendant::snp[@name='N:D377Y']/snp_call/@percent &gt; 80 and descendant::snp[@name='N:D377Y']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='N:D377Y']/snp_call/@percent &gt; 80 and descendant::snp[@name='N:D377Y']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_L18F"><xsl:choose><xsl:when test="descendant::snp[@name='S:L18F']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:L18F']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:L18F']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:L18F']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_K417T"><xsl:choose><xsl:when test="descendant::snp[@name='S:K417T']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:K417T']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:K417T']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:K417T']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_E484K"><xsl:choose><xsl:when test="descendant::snp[@name='S:E484K']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:E484K']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:E484K']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:E484K']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_N501Y"><xsl:choose><xsl:when test="descendant::snp[@name='S:N501Y']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:N501Y']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:N501Y']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:N501Y']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_H655Y"><xsl:choose><xsl:when test="descendant::snp[@name='S:H655Y']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:H655Y']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:H655Y']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:H655Y']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_A570D"><xsl:choose><xsl:when test="descendant::snp[@name='S:A570D']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:A570D']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:A570D']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:A570D']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_S982A"><xsl:choose><xsl:when test="descendant::snp[@name='S:S982A']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:S982A']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:S982A']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:S982A']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_D1118H"><xsl:choose><xsl:when test="descendant::snp[@name='S:D1118H']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:D1118H']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:D1118H']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:D1118H']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_K417N"><xsl:choose><xsl:when test="descendant::snp[@name='S:K417N']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:K417N']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:K417N']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:K417N']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_A701V"><xsl:choose><xsl:when test="descendant::snp[@name='S:A701V']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:A701V']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:A701V']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:A701V']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_K378E"><xsl:choose><xsl:when test="descendant::snp[@name='S:K378E']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:K378E']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:K378E']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:K378E']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_K378N"><xsl:choose><xsl:when test="descendant::snp[@name='S:K378N']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:K378N']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:K378N']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:K378N']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_G446D"><xsl:choose><xsl:when test="descendant::snp[@name='S:G446D']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:G446D']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:G446D']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:G446D']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_Q498R"><xsl:choose><xsl:when test="descendant::snp[@name='S:Q498R']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:Q498R']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:Q498R']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:Q498R']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_D614G"><xsl:choose><xsl:when test="descendant::snp[@name='S:D614G']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:D614G']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:D614G']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:D614G']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_P1427L"><xsl:choose><xsl:when test="descendant::snp[@name='ORF1b:P1427L']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF1b:P1427L']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='ORF1b:P1427L']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF1b:P1427L']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_Y1464C"><xsl:choose><xsl:when test="descendant::snp[@name='ORF1b:Y1464C']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF1b:Y1464C']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='ORF1b:Y1464C']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF1b:Y1464C']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_I3731V"><xsl:choose><xsl:when test="descendant::snp[@name='ORF1a:I3731V']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF1a:I3731V']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='ORF1a:I3731V']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF1a:I3731V']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_A2529V"><xsl:choose><xsl:when test="descendant::snp[@name='ORF1a:A2529V']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF1a:A2529V']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='ORF1a:A2529V']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF1a:A2529V']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_T95I"><xsl:choose><xsl:when test="descendant::snp[@name='S:T95I']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:T95I']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:T95I']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:T95I']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_A222V"><xsl:choose><xsl:when test="descendant::snp[@name='S:A222V']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:A222V']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:A222V']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:A222V']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_V1104L"><xsl:choose><xsl:when test="descendant::snp[@name='S:V1104L']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:V1104L']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:V1104L']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:V1104L']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_E239Q"><xsl:choose><xsl:when test="descendant::snp[@name='ORF3a:E239Q']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF3a:E239Q']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='ORF3a:E239Q']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF3a:E239Q']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_G38syn"><xsl:choose><xsl:when test="descendant::snp[@name='ORF7a:G38syn']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF7a:G38syn']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='ORF7a:G38syn']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF7a:G38syn']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_V1264L"><xsl:choose><xsl:when test="descendant::snp[@name='S:V1264L']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:V1264L']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:V1264L']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:V1264L']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_K48N"><xsl:choose><xsl:when test="descendant::snp[@name='ORF6:K48N']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF6:K48N']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='ORF6:K48N']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF6:K48N']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_Q2635H"><xsl:choose><xsl:when test="descendant::snp[@name='ORF1b:Q2635H']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF1b:Q2635H']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='ORF1b:Q2635H']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF1b:Q2635H']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_V71I"><xsl:choose><xsl:when test="descendant::snp[@name='ORF7a:V71I']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF7a:V71I']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='ORF7a:V71I']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF7a:V71I']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_F1925syn"><xsl:choose><xsl:when test="descendant::snp[@name='ORF1a:F1925syn']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF1a:F1925syn']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='ORF1a:F1925syn']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF1a:F1925syn']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_A2554V"><xsl:choose><xsl:when test="descendant::snp[@name='ORF1a:A2554V']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF1a:A2554V']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='ORF1a:A2554V']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF1a:A2554V']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_I1257V"><xsl:choose><xsl:when test="descendant::snp[@name='ORF1b:I1257V']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF1b:I1257V']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='ORF1b:I1257V']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF1b:I1257V']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_I882syn"><xsl:choose><xsl:when test="descendant::snp[@name='S:I882syn']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:I882syn']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:I882syn']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:I882syn']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_P314L"><xsl:choose><xsl:when test="descendant::snp[@name='ORF1b:P314L']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF1b:P314L']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='ORF1b:P314L']/snp_call/@percent &gt; 80 and descendant::snp[@name='ORF1b:P314L']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="COLOR"><xsl:choose><xsl:when test="assay/amplicon/breadth &lt; 90">#888888</xsl:when><xsl:otherwise>#000000</xsl:otherwise></xsl:choose></xsl:variable>
                    <tr>
                    <td nowrap="true" class="fixed-col-sample"><a href="{/analysis/@run_name}/{./@name}.html"><xsl:value-of select="@name"/></a></td>
	    	    <td align="center" style="color:{$COLOR}" class="fixed-col-breadth"><xsl:value-of select='format-number(assay/amplicon/breadth, "##.##")'/>%</td>
                    <td style="background:{$BG_T19R}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:T19R']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:T19R']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:T19R"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_G142D}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:G142D']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:G142D']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:G142D"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_L452R}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:L452R']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:L452R']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:L452R"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_T478K}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:T478K']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:T478K']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:T478K"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_P681R}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:P681R']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:P681R']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:P681R"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_D950N}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:D950N']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:D950N']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:D950N"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_S26L}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='ORF3a:S26L']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='ORF3a:S26L']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="ORF3a:S26L"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_I82T}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='M:I82T']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='M:I82T']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="M:I82T"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_V82A}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='ORF7a:V82A']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='ORF7a:V82A']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="ORF7a:V82A"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_T120I}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='ORF7a:T120I']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='ORF7a:T120I']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="ORF7a:T120I"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_D63G}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='N:D63G']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='N:D63G']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="N:D63G"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_R203M}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='N:R203M']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='N:R203M']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="N:R203M"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_D377Y}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='N:D377Y']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='N:D377Y']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="N:D377Y"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_I3731V}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='ORF1a:I3731V']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='ORF1a:I3731V']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="ORF1a:I3731V"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_A2529V}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='ORF1a:A2529V']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='ORF1a:A2529V']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="ORF1a:A2529V"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_T95I}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:T95I']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:T95I']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:T95I"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_V1104L}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:V1104L']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:V1104L']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:V1104L"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_A222V}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:A222V']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:A222V']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:A222V"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_E239Q}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='ORF3a:E239Q']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='ORF3a:E239Q']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="ORF3a:E239Q"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_G38syn}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='ORF7a:G38syn']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='ORF7a:G38syn']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="ORF7a:G38syn"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_A222V}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:A222V']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:A222V']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:A222V"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_V1264L}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:V1264L']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:V1264L']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:V1264L"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_K48N}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='ORF6:K48N']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='ORF6:K48N']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="ORF6:K48N"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_Q2635H}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='ORF1b:Q2635H']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='ORF1b:Q2635H']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="ORF1b:Q2635H"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_V71I}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='ORF7a:V71I']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='ORF7a:V71I']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="ORF7a:V71I"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_F1925syn}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='ORF1a:F1925syn']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='ORF1a:F1925syn']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="ORF1a:F1925syn"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_A2554V}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='ORF1a:A2554V']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='ORF1a:A2554V']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="ORF1a:A2554V"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_I1257V}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='ORF1b:I1257V']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='ORF1b:I1257V']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="ORF1b:I1257V"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_I882syn}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:I882syn']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:I882syn']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:I882syn"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_L18F}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:L18F']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:L18F']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:L18F"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_K417T}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:K417T']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:K417T']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:K417T"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_E484K}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:E484K']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:E484K']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:E484K"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_N501Y}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:N501Y']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:N501Y']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:N501Y"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_H655Y}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:H655Y']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:H655Y']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:H655Y"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_N501Y}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:N501Y']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:N501Y']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:N501Y"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_A570D}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:A570D']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:A570D']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:A570D"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_S982A}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:S982A']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:S982A']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:S982A"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_D1118H}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:D1118H']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:D1118H']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:D1118H"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_K417N}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:K417N']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:K417N']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:K417N"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_E484K}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:E484K']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:E484K']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:E484K"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_N501Y}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:N501Y']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:N501Y']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:N501Y"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_A701V}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:A701V']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:A701V']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:A701V"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_K378E}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:K378E']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:K378E']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:K378E"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_K378N}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:K378N']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:K378N']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:K378N"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_G446D}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:G446D']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:G446D']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:G446D"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_Q498R}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:Q498R']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:Q498R']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:Q498R"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_D614G}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:D614G']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:D614G']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:D614G"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_P314L}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='ORF1b:P314L']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='ORF1b:P314L']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="ORF1b:P314L"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_P1427L}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='ORF1b:P1427L']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='ORF1b:P1427L']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="ORF1b:P1427L"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_Y1464C}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='ORF1b:Y1464C']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='ORF1b:Y1464C']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="ORF1b:Y1464C"]/snp_call/@percent, "##.##")'/>%)</td>
	    	    <td nowrap="true" style="color:{$COLOR}"><xsl:for-each select="assay/amplicon/snp">
                        <xsl:if test="@name='unknown' and snp_call='_'">
	    		<xsl:value-of select="@reference"/><xsl:value-of select="@position"/>del(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%) <xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/><br/>
                        </xsl:if>
	    	    </xsl:for-each></td>
                    </tr>
                </xsl:for-each>
            </table>
            <br/>
            <a href="{@run_name}.html">Click here to return to the run summary</a>
        </body>
        </html>

    </exsl:document>
    </xsl:template>
    
    <xsl:template match="sample">
	<exsl:document method="html" href="{/analysis/@run_name}/{@name}.html">
	    <html>
	    <head>
	    	<title>ASAP Results for Sample: <xsl:value-of select="@name"/></title>
	    </head>
	    <body>
	        <center><h1>ASAP Results for Sample: <xsl:value-of select="@name"/></h1></center>
	        <br />
	        Total reads: <xsl:value-of select="@mapped_reads + @unmapped_reads"/><br/>
	        Mapped reads: <xsl:value-of select="@mapped_reads"/><br/>
	        Unmapped reads: <xsl:value-of select="@unmapped_reads"/><br/>
	        <br/>
	    	<table border="2" rules="rows">
	    	    <tr><th colspan="2">Clinical Summary for Sample: <xsl:value-of select="@name"/></th></tr>
	    	    <xsl:for-each select=".//significance">
                    <xsl:if test="not(./@flag)">
	    		    <tr>
	    		        <td><xsl:text disable-output-escaping="yes"><![CDATA[&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;]]></xsl:text></td>
	    		        <td><xsl:value-of select="."/></td>
	    		    </tr>
                    </xsl:if>
	            </xsl:for-each>
	    	</table>
	    	<br />
	    	<br />
	    	<h3>Mutations present in sample: <xsl:value-of select="@name"/></h3>
	    	<table border="2" width="100%">
	    		<tr>
	    		<th>Assay Name</th>
	    		<th># of Reads</th>
	    		<th>Coverage Breadth</th>
	    		<th>Average Depth</th>
	    		<th>Tracked Mutations</th>
	    		<th>Other SNPs found(% reads containing SNP)</th>
	    		</tr>
	    		<xsl:for-each select="assay">
	    		    <xsl:if test="@type = 'SNP' or @type = 'mixed' and amplicon/@reads &gt; 0">
	    		    <tr>
	    		        <td><xsl:value-of select="@name"/></td>
	    		        <td><xsl:value-of select="amplicon/@reads"/></td>
	    		        <td><xsl:value-of select='format-number(amplicon/breadth, "##.##")'/>%</td>
	    		        <td><xsl:value-of select='format-number(amplicon/average_depth, "##.###")'/></td>
	    		        <td><xsl:for-each select="amplicon/snp">
                                    <xsl:if test="@name!='unknown'">
	    		            <xsl:value-of select="@name"/> - <xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%)<br/>
                                    </xsl:if>
	    		        </xsl:for-each></td>
	    		        <!--<td><xsl:value-of select="amplicon/significance"/><xsl:if test="amplicon/significance/@flag"> (<xsl:value-of select="amplicon/significance/@flag"/>)</xsl:if></td>-->
	    		        <td><xsl:for-each select="amplicon/snp">
                                    <xsl:if test="@name='unknown'">
	    		            <xsl:value-of select="@position"/><xsl:value-of select="@reference"/>-><xsl:value-of select="snp_call"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%) <xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/><br/>
                                    </xsl:if>
	    		        </xsl:for-each></td>
	    		    </tr>
	    		    </xsl:if>
	    		</xsl:for-each>
	    	</table>
	    	<br />
	    </body>
	    </html>
	</exsl:document>
<xsl:variable name="SAMPLE">
<xsl:choose><xsl:when test="contains(@name, 'WMTS')">
<xsl:value-of select="substring(@name,1,25)"/>
</xsl:when><xsl:when test="contains(@name, 'ARTIC')">
<xsl:value-of select="substring(@name,1,26)"/>
</xsl:when><xsl:when test="contains(@name, 'Tiled-NEC') or contains(@name, 'Tiled_NEC')">
<xsl:value-of select="substring(@name,1,33)"/>
</xsl:when><xsl:when test="contains(@name, 'Tiled-TG1') or contains(@name, 'Tiled_TG1')">
<xsl:value-of select="substring(@name,1,27)"/>
</xsl:when><xsl:when test="contains(@name, 'Tiled-') or contains(@name, 'Tiled_')">
<xsl:value-of select="substring(@name,1,26)"/>
</xsl:when><xsl:when test="contains(@name, 'Tiled400')">
<xsl:value-of select="substring(@name,1,29)"/>
</xsl:when><xsl:when test="contains(@name, 'Tiled1000')">
<xsl:value-of select="substring(@name,1,30)"/>
</xsl:when><xsl:otherwise>
<xsl:value-of select="@name"/>
</xsl:otherwise></xsl:choose>
</xsl:variable>
  <xsl:if test="assay[1]/amplicon[1]/breadth &gt;= 90 and assay[1]/amplicon[1]/average_depth &gt;= 30">
<exsl:document method="text" href="{/analysis/@run_name}/{$SAMPLE}_gapfilled.fasta">
<xsl:for-each select="assay">
<xsl:text/>><xsl:value-of select="$SAMPLE"/>
<xsl:text>&#xa;</xsl:text>
<xsl:value-of select="./amplicon[1]/gapfilled_consensus_sequence"/>
<xsl:text>&#xa;</xsl:text>
</xsl:for-each>
<xsl:text/>
</exsl:document>
  </xsl:if>
  <xsl:if test="assay[1]/amplicon[1]/breadth &gt;= 80 and assay[1]/amplicon[1]/breadth &lt; 90 and assay[1]/amplicon[1]/average_depth &gt;= 30">
<exsl:document method="text" href="{/analysis/@run_name}/{$SAMPLE}_incomplete.fasta">
<xsl:for-each select="assay">
<xsl:text/>><xsl:value-of select="$SAMPLE"/>
<xsl:text>&#xa;</xsl:text>
<xsl:value-of select="./amplicon[1]/gapfilled_consensus_sequence"/>
<xsl:text>&#xa;</xsl:text>
</xsl:for-each>
<xsl:text/>
</exsl:document>
  </xsl:if>
  <xsl:if test="assay[1]/amplicon[1]/breadth &lt; 80 or assay[1]/amplicon[1]/average_depth &lt; 30">
<xsl:variable name="BREADTH" select='format-number(assay/amplicon/breadth, "0")'/>
<exsl:document method="text" href="{/analysis/@run_name}/{$SAMPLE}_{$BREADTH}percent.fasta">
<xsl:for-each select="assay">
<xsl:text/>><xsl:value-of select="$SAMPLE"/>
<xsl:text>&#xa;</xsl:text>
<xsl:value-of select="./amplicon[1]/gapfilled_consensus_sequence"/>
<xsl:text>&#xa;</xsl:text>
</xsl:for-each>
<xsl:text/>
</exsl:document>
  </xsl:if>
    </xsl:template>
</xsl:stylesheet>
