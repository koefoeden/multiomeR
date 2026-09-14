import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location("watch", Path(__file__).resolve().parents[2] / "scripts/local_validation/watch.py")
watch = importlib.util.module_from_spec(spec)
spec.loader.exec_module(watch)


class WatchTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.record = dict(run=self.directory.name, sha="a" * 40, profile="private", profile_sha="b" * 40, job="42")

    def result(self, **changes):
        value = dict(sha=self.record['sha'], profile_sha=self.record['profile_sha'], exit_code=0, validated=True)
        value.update(changes)
        Path(self.directory.name, 'result.json').write_text(json.dumps(value))

    def test_success_requires_output_checks(self):
        self.result(validated=False)
        self.assertEqual(watch.reconcile(self.record)[0], 'failure')
        self.result()
        self.assertEqual(watch.reconcile(self.record)[0], 'success')

    def test_code_and_profile_provenance(self):
        for change in (dict(sha='c'*40), dict(profile_sha='c'*40)):
            self.result(**change)
            self.assertEqual(watch.reconcile(self.record)[0], 'error')

    def test_nonzero_exit_cannot_pass(self):
        self.result(exit_code=143)
        self.assertEqual(watch.reconcile(self.record)[0], 'failure')

    @patch.object(watch, 'command')
    def test_cancelled_job_without_result(self, command):
        command.return_value = '42|CANCELLED by 1000|\n42.batch|CANCELLED|'
        self.assertEqual(watch.reconcile(self.record)[0], 'error')

    @patch.object(watch, 'command')
    def test_missing_accounting_is_not_success(self, command):
        command.return_value = ''
        self.assertEqual(watch.reconcile(self.record)[0], 'pending')
        command.return_value = '42|RUNNING|'
        self.assertEqual(watch.reconcile(self.record)[0], 'pending')

    def test_interrupted_submission_is_not_repeated(self):
        del self.record['job']
        self.assertEqual(watch.reconcile(self.record)[0], 'error')

    @patch.object(watch, 'command')
    def test_status_does_not_publish_private_details(self, command):
        self.record['run'] = '/private/patient-identifiers'
        watch.publish(self.record, 'success', 'Pipeline execution and output checks passed')
        payload = str(command.call_args_list)
        self.assertNotIn(self.record['run'], payload)
        self.assertNotIn(self.record['profile_sha'], payload)
        self.assertIn('local-validation/differential-analysis', payload)

    def test_rc_namespace(self):
        self.assertIsNotNone(watch.RC.fullmatch('v0.6.0-rc.1'))
        for tag in ('v0.6.0', 'ci-data-v2', 'v0.6.0-rc.1-extra'):
            self.assertIsNone(watch.RC.fullmatch(tag))

    @patch.object(watch, 'publish')
    @patch.object(watch, 'reconcile', return_value=('pending', 'Waiting'))
    @patch.object(watch, 'command')
    @patch.object(watch, 'submit')
    def test_polling_submits_each_profile_once(self, submit, command, reconcile, publish):
        root = Path(self.directory.name)
        (root / 'repository.git').mkdir()
        sha = 'a' * 40

        def respond(*args):
            if 'tag' in args:
                return 'v0.6.0-rc.1'
            if 'rev-parse' in args:
                return sha
            return ''

        def record_submission(root, mirror, sha, tag, profile, profile_repository, records):
            records[f'{sha}/{profile}'] = dict(sha=sha, tag=tag, profile=profile, state='pending')
            watch.save(root / 'state.json', records)

        command.side_effect = respond
        submit.side_effect = record_submission
        with patch('sys.argv', ['watch.py', '--root', str(root), '--profile-repository', str(root)]):
            watch.main()
            watch.main()
        self.assertEqual(submit.call_count, 2)

    @patch.object(watch, 'command')
    @patch.object(watch, 'submit')
    def test_moved_rc_tag_is_rejected(self, submit, command):
        root = Path(self.directory.name)
        (root / 'repository.git').mkdir()
        watch.save(root / 'state.json', {'previous/demo': dict(sha='b'*40, tag='v0.6.0-rc.1', state='success')})
        command.side_effect = lambda *args: 'v0.6.0-rc.1' if 'tag' in args else ('a'*40 if 'rev-parse' in args else '')
        with patch('sys.argv', ['watch.py', '--root', str(root), '--profile-repository', str(root)]):
            with self.assertRaisesRegex(RuntimeError, 'tag moved'):
                watch.main()
        submit.assert_not_called()


if __name__ == '__main__':
    unittest.main()
