/**
 * Data Downloads page — manage vendor sessions, kick off downloads,
 * watch progress.
 *
 * MVP: polls /jobs every 5s. WebSocket realtime push is wired up in
 * the backend (realtime:task:<job_id>); a follow-up can subscribe via
 * the existing useRealtimeConnection hook to remove the polling.
 */
import { useMemo, useState } from 'react'
import {
  Button,
  Card,
  Form,
  Input,
  Modal,
  Progress,
  Select,
  Space,
  Table,
  Tag,
  Tooltip,
  Typography,
  message,
  Switch,
} from 'antd'
import { DeleteOutlined, ReloadOutlined, KeyOutlined, PlayCircleOutlined } from '@ant-design/icons'
import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query'
import { dataDownloadService, vendorCredentialService } from '@/services/data-download.service'
import type { DownloadJob, DownloadJobCreateRequest, DownloadStatus, SessionLoginRequest } from '@/types/data-download'

const { Title, Text, Paragraph } = Typography

const STATUS_COLOR: Record<DownloadStatus, string> = {
  pending: 'default',
  running: 'processing',
  completed: 'success',
  failed: 'error',
  cancelled: 'warning',
}

export function DataDownloadsPage() {
  const qc = useQueryClient()
  const [vendor, setVendor] = useState<string>('lc_bio')
  const [credentialModalOpen, setCredentialModalOpen] = useState(false)

  // Catalog of vendors the backend has adapters for.
  const vendorsQuery = useQuery({
    queryKey: ['data-downloads', 'vendors'],
    queryFn: () => dataDownloadService.listVendors(),
  })

  // Saved credentials for the chosen vendor.
  const credsQuery = useQuery({
    queryKey: ['vendor-credentials', vendor],
    queryFn: () => vendorCredentialService.list(vendor).then((r) => r.items),
  })

  // Whether the vendor's daemon is currently up.
  const sessionQuery = useQuery({
    queryKey: ['data-downloads', 'session', vendor],
    queryFn: () => dataDownloadService.getSession(vendor),
    refetchInterval: 10_000,
  })

  // User's download jobs (poll every 5s while any are running).
  const jobsQuery = useQuery({
    queryKey: ['data-downloads', 'jobs'],
    queryFn: () => dataDownloadService.listJobs().then((r) => r.items),
    refetchInterval: (query) => {
      const items = (query.state.data ?? []) as DownloadJob[]
      return items.some((j) => j.status === 'running' || j.status === 'pending') ? 5_000 : false
    },
  })

  const openSessionMutation = useMutation({
    mutationFn: (req: SessionLoginRequest) => dataDownloadService.openSession(req),
    onSuccess: () => {
      message.success('Session opened')
      qc.invalidateQueries({ queryKey: ['data-downloads', 'session', vendor] })
    },
    onError: (e: Error & { response?: { data?: { detail?: string } } }) =>
      message.error(e.response?.data?.detail ?? e.message),
  })

  const closeSessionMutation = useMutation({
    mutationFn: () => dataDownloadService.closeSession(vendor),
    onSuccess: () => {
      message.success('Session closed')
      qc.invalidateQueries({ queryKey: ['data-downloads', 'session', vendor] })
    },
  })

  const createJobMutation = useMutation({
    mutationFn: (req: DownloadJobCreateRequest) => dataDownloadService.createJob(req),
    onSuccess: () => {
      message.success('Download started')
      qc.invalidateQueries({ queryKey: ['data-downloads', 'jobs'] })
    },
    onError: (e: Error & { response?: { data?: { detail?: string } } }) =>
      message.error(e.response?.data?.detail ?? e.message),
  })

  const cancelJobMutation = useMutation({
    mutationFn: (id: string) => dataDownloadService.cancelJob(id),
    onSuccess: () => qc.invalidateQueries({ queryKey: ['data-downloads', 'jobs'] }),
  })

  const sessionActive = sessionQuery.data?.active === true

  return (
    <div style={{ padding: 24 }}>
      <Title level={3}>Data Downloads</Title>
      <Paragraph type="secondary">
        Pull sequencing deliveries from vendors (联川 / 诺禾致源 / ...) directly into NGSmodule. Copy the obs path from
        the vendor&apos;s web download button — anything in quotes after <Text code>download</Text> — and paste it
        below.
      </Paragraph>

      <Card
        title={
          <Space>
            <span>Vendor</span>
            <Select
              size="small"
              style={{ minWidth: 160 }}
              loading={vendorsQuery.isLoading}
              value={vendor}
              onChange={setVendor}
              options={(vendorsQuery.data ?? []).map((v: string) => ({ value: v, label: v }))}
            />
            <Tag color={sessionActive ? 'success' : 'default'}>session {sessionActive ? 'active' : 'inactive'}</Tag>
          </Space>
        }
        extra={
          <Space>
            <Button icon={<KeyOutlined />} onClick={() => setCredentialModalOpen(true)}>
              Saved logins ({credsQuery.data?.length ?? 0})
            </Button>
            {sessionActive && (
              <Button danger onClick={() => closeSessionMutation.mutate()}>
                Close session
              </Button>
            )}
          </Space>
        }
        style={{ marginBottom: 16 }}
      >
        <SessionForm
          vendor={vendor}
          credentials={credsQuery.data ?? []}
          onSubmit={(req) => openSessionMutation.mutate(req)}
          submitting={openSessionMutation.isPending}
          disabled={sessionActive}
        />
      </Card>

      <Card title="New download" style={{ marginBottom: 16 }}>
        <DownloadForm
          vendor={vendor}
          disabled={!sessionActive}
          submitting={createJobMutation.isPending}
          onSubmit={(req) => createJobMutation.mutate(req)}
        />
      </Card>

      <Card
        title={`Jobs (${jobsQuery.data?.length ?? 0})`}
        extra={
          <Button
            icon={<ReloadOutlined />}
            onClick={() => qc.invalidateQueries({ queryKey: ['data-downloads', 'jobs'] })}
          >
            Refresh
          </Button>
        }
      >
        <JobsTable
          jobs={jobsQuery.data ?? []}
          loading={jobsQuery.isLoading}
          onCancel={(id) => cancelJobMutation.mutate(id)}
        />
      </Card>

      <CredentialsModal vendor={vendor} open={credentialModalOpen} onClose={() => setCredentialModalOpen(false)} />
    </div>
  )
}

// -------- session form --------

function SessionForm({
  vendor,
  credentials,
  onSubmit,
  submitting,
  disabled,
}: {
  vendor: string
  credentials: { id: string; label: string; email_preview: string }[]
  onSubmit: (req: SessionLoginRequest) => void
  submitting: boolean
  disabled: boolean
}) {
  const [credentialId, setCredentialId] = useState<string | undefined>(undefined)
  const [form] = Form.useForm<{ email: string; password: string }>()

  return (
    <Form
      layout="inline"
      form={form}
      disabled={disabled}
      onFinish={(values) => {
        if (credentialId) {
          onSubmit({ vendor, credential_id: credentialId })
        } else {
          onSubmit({ vendor, ...values })
        }
      }}
    >
      <Form.Item label="Use saved">
        <Select
          allowClear
          style={{ minWidth: 200 }}
          placeholder="(none — type below)"
          value={credentialId}
          onChange={setCredentialId}
          options={credentials.map((c) => ({
            value: c.id,
            label: `${c.label} — ${c.email_preview}`,
          }))}
        />
      </Form.Item>
      {!credentialId && (
        <>
          <Form.Item name="email" rules={[{ required: !credentialId }]}>
            <Input placeholder="email" autoComplete="off" style={{ width: 220 }} />
          </Form.Item>
          <Form.Item name="password" rules={[{ required: !credentialId }]}>
            <Input.Password placeholder="password" autoComplete="off" style={{ width: 200 }} />
          </Form.Item>
        </>
      )}
      <Form.Item>
        <Button type="primary" htmlType="submit" icon={<PlayCircleOutlined />} loading={submitting}>
          Open session
        </Button>
      </Form.Item>
    </Form>
  )
}

// -------- new download form --------

function DownloadForm({
  vendor,
  disabled,
  submitting,
  onSubmit,
}: {
  vendor: string
  disabled: boolean
  submitting: boolean
  onSubmit: (req: DownloadJobCreateRequest) => void
}) {
  const [form] = Form.useForm<{
    source_path: string
    dest_path: string
    auto_register: boolean
    project_name?: string
  }>()
  return (
    <Form
      form={form}
      layout="vertical"
      disabled={disabled}
      onFinish={(values) =>
        onSubmit({
          vendor,
          source_path: values.source_path,
          dest_path: values.dest_path,
          auto_register: values.auto_register,
          project_name: values.project_name,
        })
      }
      initialValues={{ auto_register: false }}
    >
      <Form.Item
        name="source_path"
        label="Source path (obs path from vendor's web UI)"
        rules={[{ required: true, message: 'paste the obs path' }]}
      >
        <Input placeholder="LC-P20260327047-.../2026-04-21_17:24:55/Data.tar" />
      </Form.Item>
      <Form.Item
        name="dest_path"
        label="Local destination directory"
        rules={[{ required: true, message: 'where to save' }]}
      >
        <Input placeholder="/data/lab/.../rawdata" />
      </Form.Item>
      <Form.Item name="auto_register" label="Auto-register as NGSmodule project on completion" valuePropName="checked">
        <Switch />
      </Form.Item>
      <Form.Item name="project_name" label="Project name (optional, auto-derived if blank)">
        <Input maxLength={100} />
      </Form.Item>
      <Form.Item>
        <Button type="primary" htmlType="submit" loading={submitting} disabled={disabled}>
          Start download
        </Button>
      </Form.Item>
    </Form>
  )
}

// -------- jobs table --------

function JobsTable({
  jobs,
  loading,
  onCancel,
}: {
  jobs: DownloadJob[]
  loading: boolean
  onCancel: (id: string) => void
}) {
  const columns = useMemo(
    () => [
      {
        title: 'Source',
        dataIndex: 'source_path',
        key: 'source_path',
        ellipsis: true,
        render: (s: string) => (
          <Tooltip title={s}>
            <Text code>{s.split('/').pop()}</Text>
          </Tooltip>
        ),
      },
      {
        title: 'Destination',
        dataIndex: 'dest_path',
        key: 'dest_path',
        ellipsis: true,
        width: 220,
      },
      {
        title: 'Status',
        dataIndex: 'status',
        key: 'status',
        width: 110,
        render: (s: DownloadStatus) => <Tag color={STATUS_COLOR[s]}>{s}</Tag>,
      },
      {
        title: 'Progress',
        dataIndex: 'progress_pct',
        key: 'progress_pct',
        width: 200,
        render: (pct: number, row: DownloadJob) => (
          <Progress
            percent={Math.round(pct)}
            status={row.status === 'failed' ? 'exception' : row.status === 'completed' ? 'success' : 'active'}
            size="small"
          />
        ),
      },
      {
        title: 'Started',
        dataIndex: 'started_at',
        key: 'started_at',
        width: 170,
        render: (v: string | null) => (v ? new Date(v).toLocaleString() : '—'),
      },
      {
        title: '',
        key: 'actions',
        width: 80,
        render: (_: unknown, row: DownloadJob) =>
          row.status === 'running' || row.status === 'pending' ? (
            <Button size="small" danger icon={<DeleteOutlined />} onClick={() => onCancel(row.id)}>
              Cancel
            </Button>
          ) : null,
      },
    ],
    [onCancel],
  )
  return (
    <Table
      rowKey="id"
      size="small"
      loading={loading}
      dataSource={jobs}
      columns={columns}
      pagination={{ pageSize: 10 }}
    />
  )
}

// -------- saved credentials modal --------

function CredentialsModal({ vendor, open, onClose }: { vendor: string; open: boolean; onClose: () => void }) {
  const qc = useQueryClient()
  const [form] = Form.useForm<{ email: string; password: string; label: string }>()

  const credsQuery = useQuery({
    queryKey: ['vendor-credentials', vendor],
    queryFn: () => vendorCredentialService.list(vendor).then((r) => r.items),
    enabled: open,
  })

  const createMutation = useMutation({
    mutationFn: (req: { vendor: string; email: string; password: string; label: string }) =>
      vendorCredentialService.create(req),
    onSuccess: () => {
      message.success('Saved')
      form.resetFields()
      qc.invalidateQueries({ queryKey: ['vendor-credentials', vendor] })
    },
  })

  const deleteMutation = useMutation({
    mutationFn: (id: string) => vendorCredentialService.delete(id),
    onSuccess: () => qc.invalidateQueries({ queryKey: ['vendor-credentials', vendor] }),
  })

  return (
    <Modal open={open} onCancel={onClose} title={`Saved logins (${vendor})`} footer={null} width={620}>
      <Table
        rowKey="id"
        size="small"
        dataSource={credsQuery.data ?? []}
        pagination={false}
        columns={[
          { title: 'Label', dataIndex: 'label' },
          { title: 'Email', dataIndex: 'email_preview' },
          {
            title: '',
            width: 70,
            render: (_, row) => (
              <Button size="small" danger icon={<DeleteOutlined />} onClick={() => deleteMutation.mutate(row.id)} />
            ),
          },
        ]}
        style={{ marginBottom: 16 }}
      />
      <Form
        form={form}
        layout="inline"
        onFinish={(v) => createMutation.mutate({ vendor, ...v, label: v.label || 'default' })}
      >
        <Form.Item name="label">
          <Input placeholder="label (default)" style={{ width: 130 }} />
        </Form.Item>
        <Form.Item name="email" rules={[{ required: true }]}>
          <Input placeholder="email" style={{ width: 180 }} autoComplete="off" />
        </Form.Item>
        <Form.Item name="password" rules={[{ required: true }]}>
          <Input.Password placeholder="password" style={{ width: 160 }} autoComplete="off" />
        </Form.Item>
        <Form.Item>
          <Button type="primary" htmlType="submit" loading={createMutation.isPending}>
            Save
          </Button>
        </Form.Item>
      </Form>
    </Modal>
  )
}
