---
layout: page
show_meta: false
title: "Numerical Packages Agenda"
header:
   image_fullwidth: "q_center_main.jpg"
permalink: "/agenda/"
---

### [Be sure to see these VIP talks]({{site.url}}{{site.baseurl}}/vip_talks/)
<br>

### [Session Selection Survey](https://forms.gle/FnEBdnpseP9xPYkh7)
Please fill out the Session Selection Survey to help us arrange
rooms for the presentations and ensure everyone has a seat.
<br>
<br>

{% include agenda %}
<br>

{% comment %}
{% for example in site.data.pathways %}
**{{example.title}}:**

   {% for agitm in site.data.agenda %}
      {% if example.session1 == agitm.id %}
* 10:30am -- [{{agitm.title|replace:"<br>"," "}}]({{site.url}}{{site.baseurl}}/session_synopses/#{{agitm.title|downcase|replace: " ", "-"|replace: "&", ""|replace: "<br>", "-"|replace: "(", ""|replace: ")", ""|replace: "/", ""}}) -- [{{agitm.vname|split: " "|join: "<br>"}}]
      {% endif %}
   {% endfor %}
   {% for agitm in site.data.agenda %}
      {% if example.session2 == agitm.id %}
* 11:45am -- [{{agitm.title|replace:"<br>"," "}}]({{site.url}}{{site.baseurl}}/session_synopses/#{{agitm.title|downcase|replace: " ", "-"|replace: "&", ""|replace: "<br>", "-"|replace: "(", ""|replace: ")", ""|replace: "/", ""}}) -- [{{agitm.vname|split: " "|join: "<br>"}}]
      {% endif %}
   {% endfor %}
   {% for agitm in site.data.agenda %}
      {% if example.session3 == agitm.id %}
* 02:35pm -- [{{agitm.title|replace:"<br>"," "}}]({{site.url}}{{site.baseurl}}/session_synopses/#{{agitm.title|downcase|replace: " ", "-"|replace: "&", ""|replace: "<br>", "-"|replace: "(", ""|replace: ")", ""|replace: "/", ""}}) -- [{{agitm.vname|split: " "|join: "<br>"}}]
      {% endif %}
   {% endfor %}
   {% for agitm in site.data.agenda %}
      {% if example.session4 == agitm.id %}
* 03:40pm -- [{{agitm.title|replace:"<br>"," "}}]({{site.url}}{{site.baseurl}}/session_synopses/#{{agitm.title|downcase|replace: " ", "-"|replace: "&", ""|replace: "<br>", "-"|replace: "(", ""|replace: ")", ""|replace: "/", ""}}) -- [{{agitm.vname|split: " "|join: "<br>"}}]
      {% endif %}
   {% endfor %}
<br>
{% endfor %}
{% endcomment %}

### [ATPESC 2025 Main Agenda](https://extremecomputingtraining.anl.gov/agenda-2025/)

{% include link-shortcuts %}


