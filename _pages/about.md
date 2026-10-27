---
permalink: /
layout: home
title: "Computational chemist building AI-powered tools for drug discovery"
author_profile: false
redirect_from:
  - /about/
  - /about.html
---

<div class="home-hero">
  <p class="home-kicker">CHEMINFORMATICS · AI · MOLECULAR DESIGN</p>
  <h1>Computational chemistry for better drug discovery.</h1>
  <p class="home-lede">I develop cheminformatics workflows, machine-learning models, and accessible scientific software for molecular design, virtual screening, and ADMET prediction.</p>
  <div class="home-actions">
    <a class="btn btn--primary" href="{{ '/portfolio/' | relative_url }}">Explore projects</a>
    <a class="btn btn--inverse" href="{{ '/cv/' | relative_url }}">View my CV</a>
  </div>
</div>

## What I work on

My work connects chemistry, data science, and software engineering to make computational methods more useful in real research settings.

<div class="home-grid">
  <div class="home-card"><span class="home-card__icon">01</span><h3>Cheminformatics</h3><p>RDKit-based molecular representations, chemical databases, similarity searching, and reproducible data workflows.</p></div>
  <div class="home-card"><span class="home-card__icon">02</span><h3>AI for drug discovery</h3><p>QSAR, deep learning, ADMET prediction, and AI-assisted approaches to explore chemical space.</p></div>
  <div class="home-card"><span class="home-card__icon">03</span><h3>Molecular modeling</h3><p>Molecular docking, virtual screening, molecular dynamics, and structure-based drug design.</p></div>
</div>

## Featured work

<div class="home-projects">
  <article class="home-project"><p class="home-project__eyebrow">WEB APPLICATIONS</p><h3>Chemoinformatics toolkit</h3><p>A collection of browser-based tools for molecular analysis, ADME testing, drug-likeness filters, and chemical data exploration.</p><a href="https://chemoinfo.streamlit.app/">Open the toolkit →</a></article>
  <article class="home-project"><p class="home-project__eyebrow">EDUCATION &amp; OPEN SCIENCE</p><h3>Chemoinformatics tutorials</h3><p>Practical notebooks and examples for learning cheminformatics, machine learning, and molecular modeling.</p><a href="https://github.com/Aouidate/Chemoinformatics-tutos">Browse the tutorials →</a></article>
  <article class="home-project"><p class="home-project__eyebrow">RESEARCH SOFTWARE</p><h3>Open-source experiments</h3><p>Explore my GitHub work across RDKit, Streamlit, Django, docking, QSAR, and AI-assisted drug discovery.</p><a href="https://github.com/Aouidate">View GitHub →</a></article>
</div>

## Latest writing

{% assign recent_posts = site.posts | where_exp: "post", "post.date <= site.time" | sort: "date" | reverse %}
<div class="home-posts">
{% for post in recent_posts limit:3 %}
  <a class="home-post" href="{{ post.url | relative_url }}"><span>{{ post.date | date: "%B %Y" }}</span><strong>{{ post.title }}</strong></a>
{% endfor %}
</div>

## Let’s collaborate

I’m interested in research and software projects at the intersection of computational chemistry, machine learning, and drug discovery.

<a class="btn btn--primary" href="mailto:{{ site.author.email }}">Get in touch</a>
