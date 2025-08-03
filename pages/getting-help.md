---
layout: page
show_meta: false
title: "Getting help"
header:
   image_fullwidth: "tech_support8.jpg"
permalink: "/getting-help/"
---
{% assign vroom = nil %}
{% for vr in site.data.vrooms %}
  {% if vr.name == "Amphitheater" %}
    {% assign vroom = vr %}
    {% break %}
  {% endif %}
{% endfor %}

SLACK ALCF-Workshops workspace provides channels for help

* [#atpesc-2025-helpdesk](https://alcf-workshops.slack.com/archives/C096F4REXQE)
* [#atpesc-2025-track-5-numerical](https://alcf-workshops.slack.com/archives/C0977BUSQ72)
* [#atpesc-2025-track-5-numerical-breakout](https://alcf-workshops.slack.com/archives/C096F44M10W)

{% include link-shortcuts %}
