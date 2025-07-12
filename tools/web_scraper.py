import requests
from bs4 import BeautifulSoup
import argparse

def scrape_url(url):
    """
    Scrapes the content of a given URL and returns the text.
    """
    try:
        response = requests.get(url, timeout=15)
        response.raise_for_status()
        soup = BeautifulSoup(response.content, 'html.parser')
        return soup.get_text()
    except requests.exceptions.RequestException as e:
        return f"Error fetching URL: {e}"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Scrape the text content of a web page.")
    parser.add_argument("url", help="The URL to scrape.")
    args = parser.parse_args()

    content = scrape_url(args.url)
    print(content) 