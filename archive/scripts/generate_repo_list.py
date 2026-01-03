#!/usr/bin/env python3
"""Generate organized MDC file of all GitHub repositories."""
import json
from datetime import datetime
from collections import defaultdict

with open('/tmp/repos.json', 'r') as f:
    repos = json.load(f)

repos.sort(key=lambda x: x.get('updatedAt', ''), reverse=True)

forks = [r for r in repos if r.get('isFork', False)]
own_repos = [r for r in repos if not r.get('isFork', False)]

def format_date(date_str):
    if not date_str:
        return "Unknown"
    try:
        dt = datetime.fromisoformat(date_str.replace('Z', '+00:00'))
        return dt.strftime('%Y-%m-%d')
    except:
        return date_str[:10] if date_str else "Unknown"

# Group by language
by_language = defaultdict(list)
for r in repos:
    lang = r.get('primaryLanguage', {}).get('name', 'Unknown') if r.get('primaryLanguage') else 'Unknown'
    by_language[lang].append(r)

output = []

# Top languages
output.append("## 🗣️ Top Languages by Count")
output.append("")
lang_counts = sorted(by_language.items(), key=lambda x: len(x[1]), reverse=True)[:15]
for lang, lang_repos in lang_counts:
    output.append(f"- **{lang}:** {len(lang_repos)} repos")
output.append("")

# Your Repos
output.append("## ⭐ Your Repositories (Not Forks)")
output.append("")
output.append(f"**Total: {len(own_repos)} repos**")
output.append("")
for r in own_repos:
    name = r.get('name', 'Unknown')
    desc = r.get('description', 'No description') or 'No description'
    updated = format_date(r.get('updatedAt', ''))
    created = format_date(r.get('createdAt', ''))
    stars = r.get('stargazerCount', 0)
    lang = r.get('primaryLanguage', {}).get('name', 'Unknown') if r.get('primaryLanguage') else 'Unknown'
    url = r.get('url', '')
    archived = "🔒 ARCHIVED" if r.get('isArchived', False) else ""
    
    output.append(f"### {name} {archived}")
    output.append(f"- **Description:** {desc}")
    output.append(f"- **Language:** {lang} | **Stars:** {stars}")
    output.append(f"- **Created:** {created} | **Updated:** {updated}")
    output.append(f"- **URL:** {url}")
    output.append("")

# Forks summary
output.append("## 🍴 Forks")
output.append("")
output.append(f"**Total: {len(forks)} forks**")
output.append("")
output.append("### Recent Forks (Top 50)")
output.append("")
for r in forks[:50]:
    name = r.get('name', 'Unknown')
    desc = r.get('description', 'No description') or 'No description'
    updated = format_date(r.get('updatedAt', ''))
    lang = r.get('primaryLanguage', {}).get('name', 'Unknown') if r.get('primaryLanguage') else 'Unknown'
    url = r.get('url', '')
    
    output.append(f"- **{name}** ({lang}) - {desc[:70]} | Updated: {updated} | [View]({url})")

if len(forks) > 50:
    output.append("")
    output.append(f"*... and {len(forks) - 50} more forks*")

# Complete list
output.append("")
output.append("## 📋 Complete Repository List (All 405)")
output.append("")
output.append("*Sorted by last updated (most recent first)*")
output.append("")
for i, r in enumerate(repos, 1):
    name = r.get('name', 'Unknown')
    desc = r.get('description', 'No description') or 'No description'
    updated = format_date(r.get('updatedAt', ''))
    lang = r.get('primaryLanguage', {}).get('name', 'Unknown') if r.get('primaryLanguage') else 'Unknown'
    url = r.get('url', '')
    is_fork = "🍴" if r.get('isFork', False) else "⭐"
    archived = "🔒" if r.get('isArchived', False) else ""
    
    output.append(f"{i}. {is_fork} {archived} **{name}** ({lang}) - {desc[:60]} | Updated: {updated} | [View]({url})")

# Append to file
with open('.cursor/GITHUB_REPOSITORIES.md', 'a') as f:
    f.write('\n'.join(output))

print(f"✅ Generated .cursor/GITHUB_REPOSITORIES.md")
print(f"   - {len(own_repos)} own repos")
print(f"   - {len(forks)} forks")
print(f"   - {len(repos)} total repos")

"""Generate organized MDC file of all GitHub repositories."""
import json
from datetime import datetime
from collections import defaultdict

with open('/tmp/repos.json', 'r') as f:
    repos = json.load(f)

repos.sort(key=lambda x: x.get('updatedAt', ''), reverse=True)

forks = [r for r in repos if r.get('isFork', False)]
own_repos = [r for r in repos if not r.get('isFork', False)]

def format_date(date_str):
    if not date_str:
        return "Unknown"
    try:
        dt = datetime.fromisoformat(date_str.replace('Z', '+00:00'))
        return dt.strftime('%Y-%m-%d')
    except:
        return date_str[:10] if date_str else "Unknown"

# Group by language
by_language = defaultdict(list)
for r in repos:
    lang = r.get('primaryLanguage', {}).get('name', 'Unknown') if r.get('primaryLanguage') else 'Unknown'
    by_language[lang].append(r)

output = []

# Top languages
output.append("## 🗣️ Top Languages by Count")
output.append("")
lang_counts = sorted(by_language.items(), key=lambda x: len(x[1]), reverse=True)[:15]
for lang, lang_repos in lang_counts:
    output.append(f"- **{lang}:** {len(lang_repos)} repos")
output.append("")

# Your Repos
output.append("## ⭐ Your Repositories (Not Forks)")
output.append("")
output.append(f"**Total: {len(own_repos)} repos**")
output.append("")
for r in own_repos:
    name = r.get('name', 'Unknown')
    desc = r.get('description', 'No description') or 'No description'
    updated = format_date(r.get('updatedAt', ''))
    created = format_date(r.get('createdAt', ''))
    stars = r.get('stargazerCount', 0)
    lang = r.get('primaryLanguage', {}).get('name', 'Unknown') if r.get('primaryLanguage') else 'Unknown'
    url = r.get('url', '')
    archived = "🔒 ARCHIVED" if r.get('isArchived', False) else ""
    
    output.append(f"### {name} {archived}")
    output.append(f"- **Description:** {desc}")
    output.append(f"- **Language:** {lang} | **Stars:** {stars}")
    output.append(f"- **Created:** {created} | **Updated:** {updated}")
    output.append(f"- **URL:** {url}")
    output.append("")

# Forks summary
output.append("## 🍴 Forks")
output.append("")
output.append(f"**Total: {len(forks)} forks**")
output.append("")
output.append("### Recent Forks (Top 50)")
output.append("")
for r in forks[:50]:
    name = r.get('name', 'Unknown')
    desc = r.get('description', 'No description') or 'No description'
    updated = format_date(r.get('updatedAt', ''))
    lang = r.get('primaryLanguage', {}).get('name', 'Unknown') if r.get('primaryLanguage') else 'Unknown'
    url = r.get('url', '')
    
    output.append(f"- **{name}** ({lang}) - {desc[:70]} | Updated: {updated} | [View]({url})")

if len(forks) > 50:
    output.append("")
    output.append(f"*... and {len(forks) - 50} more forks*")

# Complete list
output.append("")
output.append("## 📋 Complete Repository List (All 405)")
output.append("")
output.append("*Sorted by last updated (most recent first)*")
output.append("")
for i, r in enumerate(repos, 1):
    name = r.get('name', 'Unknown')
    desc = r.get('description', 'No description') or 'No description'
    updated = format_date(r.get('updatedAt', ''))
    lang = r.get('primaryLanguage', {}).get('name', 'Unknown') if r.get('primaryLanguage') else 'Unknown'
    url = r.get('url', '')
    is_fork = "🍴" if r.get('isFork', False) else "⭐"
    archived = "🔒" if r.get('isArchived', False) else ""
    
    output.append(f"{i}. {is_fork} {archived} **{name}** ({lang}) - {desc[:60]} | Updated: {updated} | [View]({url})")

# Append to file
with open('.cursor/GITHUB_REPOSITORIES.md', 'a') as f:
    f.write('\n'.join(output))

print(f"✅ Generated .cursor/GITHUB_REPOSITORIES.md")
print(f"   - {len(own_repos)} own repos")
print(f"   - {len(forks)} forks")
print(f"   - {len(repos)} total repos")

"""Generate organized MDC file of all GitHub repositories."""
import json
from datetime import datetime
from collections import defaultdict

with open('/tmp/repos.json', 'r') as f:
    repos = json.load(f)

repos.sort(key=lambda x: x.get('updatedAt', ''), reverse=True)

forks = [r for r in repos if r.get('isFork', False)]
own_repos = [r for r in repos if not r.get('isFork', False)]

def format_date(date_str):
    if not date_str:
        return "Unknown"
    try:
        dt = datetime.fromisoformat(date_str.replace('Z', '+00:00'))
        return dt.strftime('%Y-%m-%d')
    except:
        return date_str[:10] if date_str else "Unknown"

# Group by language
by_language = defaultdict(list)
for r in repos:
    lang = r.get('primaryLanguage', {}).get('name', 'Unknown') if r.get('primaryLanguage') else 'Unknown'
    by_language[lang].append(r)

output = []

# Top languages
output.append("## 🗣️ Top Languages by Count")
output.append("")
lang_counts = sorted(by_language.items(), key=lambda x: len(x[1]), reverse=True)[:15]
for lang, lang_repos in lang_counts:
    output.append(f"- **{lang}:** {len(lang_repos)} repos")
output.append("")

# Your Repos
output.append("## ⭐ Your Repositories (Not Forks)")
output.append("")
output.append(f"**Total: {len(own_repos)} repos**")
output.append("")
for r in own_repos:
    name = r.get('name', 'Unknown')
    desc = r.get('description', 'No description') or 'No description'
    updated = format_date(r.get('updatedAt', ''))
    created = format_date(r.get('createdAt', ''))
    stars = r.get('stargazerCount', 0)
    lang = r.get('primaryLanguage', {}).get('name', 'Unknown') if r.get('primaryLanguage') else 'Unknown'
    url = r.get('url', '')
    archived = "🔒 ARCHIVED" if r.get('isArchived', False) else ""
    
    output.append(f"### {name} {archived}")
    output.append(f"- **Description:** {desc}")
    output.append(f"- **Language:** {lang} | **Stars:** {stars}")
    output.append(f"- **Created:** {created} | **Updated:** {updated}")
    output.append(f"- **URL:** {url}")
    output.append("")

# Forks summary
output.append("## 🍴 Forks")
output.append("")
output.append(f"**Total: {len(forks)} forks**")
output.append("")
output.append("### Recent Forks (Top 50)")
output.append("")
for r in forks[:50]:
    name = r.get('name', 'Unknown')
    desc = r.get('description', 'No description') or 'No description'
    updated = format_date(r.get('updatedAt', ''))
    lang = r.get('primaryLanguage', {}).get('name', 'Unknown') if r.get('primaryLanguage') else 'Unknown'
    url = r.get('url', '')
    
    output.append(f"- **{name}** ({lang}) - {desc[:70]} | Updated: {updated} | [View]({url})")

if len(forks) > 50:
    output.append("")
    output.append(f"*... and {len(forks) - 50} more forks*")

# Complete list
output.append("")
output.append("## 📋 Complete Repository List (All 405)")
output.append("")
output.append("*Sorted by last updated (most recent first)*")
output.append("")
for i, r in enumerate(repos, 1):
    name = r.get('name', 'Unknown')
    desc = r.get('description', 'No description') or 'No description'
    updated = format_date(r.get('updatedAt', ''))
    lang = r.get('primaryLanguage', {}).get('name', 'Unknown') if r.get('primaryLanguage') else 'Unknown'
    url = r.get('url', '')
    is_fork = "🍴" if r.get('isFork', False) else "⭐"
    archived = "🔒" if r.get('isArchived', False) else ""
    
    output.append(f"{i}. {is_fork} {archived} **{name}** ({lang}) - {desc[:60]} | Updated: {updated} | [View]({url})")

# Append to file
with open('.cursor/GITHUB_REPOSITORIES.md', 'a') as f:
    f.write('\n'.join(output))

print(f"✅ Generated .cursor/GITHUB_REPOSITORIES.md")
print(f"   - {len(own_repos)} own repos")
print(f"   - {len(forks)} forks")
print(f"   - {len(repos)} total repos")



