
SCRAPING PUBMED TO ANALYZE ENCODE PROJECT PUBLICATIONS

Introduction:

The Encyclopedia Of DNA Elements (ENCODE) is a colaboration of research groups funded by the National Human Genome Research Institute(NHGRI) in 2003. It was planned as a follow-up to the Human Genome Project and has already begun its fourth phase on February 2017.

The aim of this project is to identify all functional elements in the human genome sequence (only the role of 1.5% of the genome is known) and build a list of functional elements in the human genome. The ENCODE Project has given researchers insight into how the human genome functions.
All the abovementioned serve as an expansive resource for understanding human health and disease.

A very important point for researchers is to be able to access the information from databases like PubMed, GEO, etc. Entrez is a data retrieval system that provides users access to NCBI’s databases manually, or by using Biopython’s Bio.Entrez module wich allows you to search PubMed or download GenBank records from within a Python script.


Objective:

To assess:
The number of ENCODE publications per year, the number of ENCODE publications over time and the journal most published in

To create several barplots to show our results

Methods:

In this project, we use the Bio.Entrez module to to obtain journal information from PubMed on the basis of a csv file where we had the PMIDs numbers from 690 ENTREZ papers. Our variables of interest were journal name, year of publication and name of the last author.

For the barplot creation we imported the pandas library, numpy, plotly.graph_objs, plotly.plotly, and matplotlib.pyplot


Code:


```python
#Import modules
import pandas as pd
import numpy as np
from Bio import Entrez
import matplotlib.pyplot as plt
import plotly
import plotly.plotly as py
import plotly.graph_objs as go
from collections import Counter
%matplotlib inline 

#Register with Entrez
Entrez.email = "maria.casal.dominguez@gmail.com"

#Define a function to grab the desired attributes of the PMIDs from PubMed
def get_pubmed_data(pmids, attrb_list):
   
    """Returns full PubMed data records for desired PMIDs in XML format. PMIDs can be found online in PubMed
    and can be accepted individually or as a list. Desired data from PMIDs ('attrb_list')
    can be viewed here *insert text file's github link here* """
        
    pubs_list = []
    pmid_number = len(pmids.index)
    iteration_number = 1

    while True:
        if (iteration_number*200 > pmid_number):
            upper_limit = pmid_number
        else:
            upper_limit = iteration_number*200

        #Convert the CSV file into a list
        pmids_list = pmids["PMID"][iteration_number*200-200:upper_limit].to_csv(index=False)

        xml_doc = Entrez.efetch(db='pubmed', id=pmids_list, retmode='xml', rettype='docsum')

        for pub in Entrez.parse(xml_doc):
            pub_attribs = []
            for attrib in attrb_list:
                pub_attribs.append(pub[attrib])
            pubs_list.append(pub_attribs)

        if (iteration_number*200 > pmid_number):
            break
        else:
            iteration_number += 1

    pd_docs = pd.DataFrame(pubs_list, columns=attrb_list)
    
    #Clean the data
    pd_docs.columns = ['PMID', 'Full Journal Name', 'Last Author', 'Publication Date']
    pd_docs['Publication Date'] = pd_docs['Publication Date'].apply(lambda x: str(x)[:4])
    pd_docs['Full Journal Name'] = pd_docs['Full Journal Name'].apply(lambda x: str(x)[:40])
    pd_docs["Author: Journal"] = pd_docs["Last Author"] + ": " + pd_docs["Full Journal Name"]

    #Show the DataFrame
    display(pd_docs.head(10))
    print(pd_docs.shape)

    #Export to .csv
    pd_docs.to_csv("pd_docs.csv")
    
    return pd_docs
```


```python
#Import the csv with the pmids as a DataFrame
pmids = pd.read_csv("publications.csv")
  
#List desired PMID attributes to get from PubMed
attrb_list = ["Id", "FullJournalName", "LastAuthor", "PubDate"]

#Call the function to grab the data
pd_docs = get_pubmed_data(pmids, attrb_list)
```


<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>PMID</th>
      <th>Full Journal Name</th>
      <th>Last Author</th>
      <th>Publication Date</th>
      <th>Author: Journal</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>18665130</td>
      <td>Nature genetics</td>
      <td>Gingeras TR</td>
      <td>2008</td>
      <td>Gingeras TR: Nature genetics</td>
    </tr>
    <tr>
      <th>1</th>
      <td>17568007</td>
      <td>Genome research</td>
      <td>Stamatoyannopoulos JA</td>
      <td>2007</td>
      <td>Stamatoyannopoulos JA: Genome research</td>
    </tr>
    <tr>
      <th>2</th>
      <td>17166863</td>
      <td>Nucleic acids research</td>
      <td>Kent WJ</td>
      <td>2007</td>
      <td>Kent WJ: Nucleic acids research</td>
    </tr>
    <tr>
      <th>3</th>
      <td>17567993</td>
      <td>Genome research</td>
      <td>Gerstein MB</td>
      <td>2007</td>
      <td>Gerstein MB: Genome research</td>
    </tr>
    <tr>
      <th>4</th>
      <td>17567995</td>
      <td>Genome research</td>
      <td>Sidow A</td>
      <td>2007</td>
      <td>Sidow A: Genome research</td>
    </tr>
    <tr>
      <th>5</th>
      <td>18258921</td>
      <td>Genome research</td>
      <td>Liu XS</td>
      <td>2008</td>
      <td>Liu XS: Genome research</td>
    </tr>
    <tr>
      <th>6</th>
      <td>17568011</td>
      <td>Genome research</td>
      <td>Baxevanis AD</td>
      <td>2007</td>
      <td>Baxevanis AD: Genome research</td>
    </tr>
    <tr>
      <th>7</th>
      <td>19425134</td>
      <td>Genome informatics. International Confer</td>
      <td>Tullius TD</td>
      <td>2008</td>
      <td>Tullius TD: Genome informatics. International ...</td>
    </tr>
    <tr>
      <th>8</th>
      <td>21439813</td>
      <td>Current opinion in structural biology</td>
      <td>Tullius TD</td>
      <td>2011</td>
      <td>Tullius TD: Current opinion in structural biology</td>
    </tr>
    <tr>
      <th>9</th>
      <td>19286520</td>
      <td>Science (New York, N.Y.)</td>
      <td>Margulies EH</td>
      <td>2009</td>
      <td>Margulies EH: Science (New York, N.Y.)</td>
    </tr>
  </tbody>
</table>
</div>


    (691, 5)



```python
#Define a function to sort desired data
def get_variable(dataFrame, variable):
    
    """Gets desired column data from DataFrame and sorts from most to least frequent. 
    Data are assigned to 'variable_results', data frequency is assigned to 'counts'."""
      
    #Create an empty list
    frequency_list = []
   
    #Fill the list with the frequency of the data
    for i in dataFrame[variable]:
        frequency_list.append(i)
         
    #Use Counter to count how many times a journal appears
    journal_count = Counter(frequency_list)

    #Sort the frequency data by most to least common
    sorted_journal_count = journal_count.most_common()

    #Create 2 new lists for the variables that you want
    global variable_results
    variable_results = []
    global counts
    counts = []
    
    #Fill variable_results and counts
    for i in sorted_journal_count:
        key, count = i
        variable_results.append(key)
        counts.append(count)

    return variable_results, counts
```


```python
#Authenticate plotly
plotly.tools.set_credentials_file(username='mpagan2', api_key='9oJfHnTtef1NWTokn4lI')

def hbargraph(xaxis, yaxis, title):

    """Using plotly, develops a horizontal bargraph with hover-over value labels and a custom title"""
    
    data = [go.Bar(
                x= xaxis,
                y= yaxis,
                orientation = 'h',
                marker = dict(
                    color = 'rgba(0, 158, 28, 0.6)',
                    line = dict(
                        color = 'rgba(0, 158, 28, 1.0)',
                    width = 3))
                )
           ]

    layout = go.Layout(
            title = title,
            autosize=False,
            width=1000,
            height=500,
            margin=go.Margin(
                l=300,
                r=50,
                b=100,
                t=100,
                pad=4
            ),
        )

    fig = go.Figure(data=data, layout=layout)
    graph = py.iplot(fig, filename=title)
    
    return graph

#Call ‘get_variable’ to obtain the desired data to fill the horizontal bar graph with
get_variable(pd_docs, "Full Journal Name")

#Call ‘hbargraph’ to create the plot
hbargraph(counts[:10], variable_results[:10], "Top 10 Journals ENCODE Authors Publish In")

```




<iframe id="igraph" scrolling="no" style="border:none;" seamless="seamless" src="https://plot.ly/~mpagan2/8.embed" height="500px" width="1000px"></iframe>




```python
get_variable(pd_docs, "Author: Journal")
hbargraph(counts[:10], variable_results[:10], "Top 10 Journals Published In By Single ENCODE Author")
```




<iframe id="igraph" scrolling="no" style="border:none;" seamless="seamless" src="https://plot.ly/~mpagan2/12.embed" height="500px" width="1000px"></iframe>




```python
get_variable(pd_docs, "Last Author")
hbargraph(counts[:10], variable_results[:10], "Top 10 ENCODE Authors")
```




<iframe id="igraph" scrolling="no" style="border:none;" seamless="seamless" src="https://plot.ly/~mpagan2/14.embed" height="500px" width="1000px"></iframe>




```python
#Define function to make bar graph
def bargraph(xaxis, yaxis, bar_labels, title):

    """Using plotly, develops a bar graph with custom hover-over value label descriptions and a custom title"""
    
    trace0 = go.Bar(
        x = xaxis,
        y = yaxis,
        text = bar_labels,
        marker=dict(
            color='rgb(153, 153, 255)',
            line=dict(
                color='rgb(8,48,107)',
                width=1.5,
            )
        ),
        opacity=0.6
    )

    data = [trace0]
    layout = go.Layout(
        title= title,
        xaxis=dict(
            autotick=False,
            ticks='outside',
            tick0=0,
            dtick=1,
            ticklen=8,
            tickwidth=4,
            tickcolor='#000'
        ),
    )

    fig = go.Figure(data=data, layout=layout)
    graph = py.iplot(fig, filename=title)
    return graph

#Grab unique years from "PubDate" column of pd_docs, sort the years, and create "'Year' Publications" labels
get_year_labels = pd_docs["Publication Date"].unique()
sorted_year_labels = sorted(get_year_labels)
labels = [i + " Publications" for i in sorted_year_labels]

#Get x-axis years and sort chronologically
years = pd_docs["Publication Date"].unique().tolist()
sorted_years = sorted(years)

#Get y-axis values
pubs_by_year = pd_docs.groupby('Publication Date')['PMID'].nunique().tolist()

#Call the function
bargraph(sorted_years, pubs_by_year, labels, "ENCODE Publications By Year")
```

    High five! You successfully sent some data to your account on plotly. View your plot in your browser at https://plot.ly/~mpagan2/0 or inside your plot.ly account where it is named 'ENCODE Publications By Year'





<iframe id="igraph" scrolling="no" style="border:none;" seamless="seamless" src="https://plot.ly/~mpagan2/0.embed" height="525px" width="100%"></iframe>




```python
def linegraph(xaxis, yaxis, line_labels, title):

    
    """Using plotly, develops a line graph with custom hover-over value label descriptions and a custom title"""
    
    trace = go.Scatter(
        x = xaxis,
        y = yaxis,
        text = line_labels,
        marker=dict(
                color='rgb(153, 153, 255)'
        ),
    )

    line_layout = go.Layout(
            title= title,
            xaxis=dict(
                autotick=False,
                ticks='outside',
                tick0=0,
                dtick=1,
                ticklen=8,
                tickwidth=4,
                tickcolor='#000'
            ),
    )

    line = [trace]

    fig = go.Figure(data=line, layout=line_layout)
    graph = py.iplot(fig, filename='Publications Over Time')
    
    return graph

#Sum the publications sequentially by year
sum_pubs = np.cumsum(pubs_by_year)

#Grab unique years from "PubDate" column of pd_docs, sort the years, and create "Total Number of Publications in 'Year'" labels
get_total_year_labels = pd_docs["Publication Date"].unique()
sorted_total_year_labels = sorted(get_total_year_labels)
total_labels = ["Total Number of Publications in " + i for i in sorted_total_year_labels]

#Call the line graph function
linegraph(sorted_years, sum_pubs, total_labels, "ENCODE Publications Over Time")
```




<iframe id="igraph" scrolling="no" style="border:none;" seamless="seamless" src="https://plot.ly/~mpagan2/20.embed" height="525px" width="100%"></iframe>

