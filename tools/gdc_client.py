import requests
import json

GDC_API_ENDPOINT = "https://api.gdc.cancer.gov"

def query_gdc_api(endpoint, payload):
    """
    Generic function to query the GDC API.
    """
    url = f"{GDC_API_ENDPOINT}/{endpoint}"
    try:
        response = requests.post(url, json=payload)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error querying GDC API: {e}")
        return None

if __name__ == '__main__':
    # Example usage: Find the top 10 cancer types by case count
    payload = {
        "facets": "primary_site",
        "size": "10",
        "format": "json"
    }
    data = query_gdc_api("projects", payload)

    if data:
        print("Top 10 Cancer Types by Case Count (from GDC):")
        for bucket in data.get("data", {}).get("aggregations", {}).get("primary_site", {}).get("buckets", []):
            print(f"- {bucket['key']}: {bucket['doc_count']} cases") 