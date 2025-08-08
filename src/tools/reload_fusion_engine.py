import requests

def reload_fusion_engine():
    """
    Calls the /reload endpoint of the fusion-engine to trigger a client reload.
    """
    base_url = "https://crispro--fusion-engine-v1-fusionengine-api.modal.run"
    endpoint = "/reload"
    url = f"{base_url}{endpoint}"

    print(f"Sending POST request to {url} to trigger client reload...")
    
    try:
        response = requests.post(url, verify=False)
        if response.status_code == 200:
            print("Reload request successful!")
            print(response.json())
        else:
            print(f"Reload request failed with status code: {response.status_code}")
            print("Response:")
            print(response.text)
    except requests.exceptions.RequestException as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    reload_fusion_engine() 