

def results_to_text(results: dict, job_name: str, output_path: str):
    with open(f"{output_path}{job_name}.txt", 'w') as writer:
        writer.write('\n'.join([f"{k}:{v}" for k, v in results.items()]))
