import modal

medgemma_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "google-generativeai",
        "tensorflow",
        "Pillow",
        "fastapi",
        "uvicorn"
    )
) 